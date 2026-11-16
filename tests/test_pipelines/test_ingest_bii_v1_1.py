"""Offline checks for the BII v1.1 ingest: ETag recomputation, resumed streams, manifests."""

import base64
import hashlib
import io
import sys
import urllib.error
from pathlib import Path

import google_crc32c
import numpy as np
import pytest

PREPROCESSING = Path(__file__).resolve().parents[2] / "pipelines/preprocessing"
sys.path.insert(0, str(PREPROCESSING))
import ingest_bii_v1_1 as bii  # noqa: E402

MIB = bii.MIB


def s3_etag(data: bytes, part: int | None) -> str:
    if part is None:
        return hashlib.md5(data).hexdigest()
    parts = [hashlib.md5(data[i : i + part]).digest() for i in range(0, len(data), part)]
    return f"{hashlib.md5(b''.join(parts)).hexdigest()}-{len(parts)}"


def test_part_size_matches_the_real_bii_2017_etag():
    assert bii.part_size(34_763_382_112, 4145) == 8 * MIB
    with pytest.raises(ValueError):
        bii.part_size(10, 3)


@pytest.mark.parametrize("part", [None, MIB, 2 * MIB])
def test_stream_hasher_reproduces_s3_etags_across_odd_chunks(part):
    data = np.random.default_rng(0).bytes(5 * MIB + 12345)
    hasher = bii.StreamHasher(len(data), s3_etag(data, part))
    for start in range(0, len(data), 777_777):
        hasher.update(data[start : start + 777_777])
    result = hasher.verified()
    assert result["sha256"] == hashlib.sha256(data).hexdigest()
    assert result["crc32c"] == base64.b64encode(google_crc32c.Checksum(data).digest()).decode()


def test_stream_hasher_rejects_corrupted_bytes():
    data = bytes(3 * MIB)
    hasher = bii.StreamHasher(len(data), s3_etag(data, MIB))
    hasher.update(data[:-1] + b"x")
    with pytest.raises(ValueError):
        hasher.verified()


class FakeResponse(io.BytesIO):
    def __init__(self, data: bytes, etag: str, fail_after: int | None = None):
        super().__init__(data)
        self.status, self.headers, self.fail_after = 206, {"ETag": f'"{etag}"'}, fail_after

    def read(self, size: int = -1) -> bytes:
        if self.fail_after is not None and self.tell() >= self.fail_after:
            raise ConnectionResetError("dropped")
        return super().read(size)


def test_stream_source_resumes_from_the_received_offset(monkeypatch):
    data = np.random.default_rng(1).bytes(3 * MIB + 5)
    etag = s3_etag(data, MIB)
    source = {"url": "https://example/bii.tif", "size": len(data), "etag": etag}
    ranges = []

    def urlopen(request, timeout):  # noqa: ARG001
        offset = int(request.get_header("Range").split("=")[1].rstrip("-"))
        ranges.append(offset)
        if len(ranges) == 2:
            raise urllib.error.URLError("temporary")
        return FakeResponse(data[offset:], etag, fail_after=MIB if len(ranges) == 1 else None)

    monkeypatch.setattr(bii.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(bii.time, "sleep", lambda _: None)
    monkeypatch.setattr(bii, "READ_CHUNK", MIB // 2)
    received = bytearray()
    assert bii.stream_source(source, received.extend)["etag"] == etag
    assert bytes(received) == data
    assert ranges == [0, MIB, MIB]


OBJ = {
    "source_url": f"{bii.SOURCE_URL}/bii_2020.tif",
    "source_etag": "abc-2",
    "source_size": 34_763_382_112,
    "source_last_modified": "Mon, 13 Jul 2026 15:00:41 GMT",
    "gcs_uri": f"gs://{bii.BUCKET}/{bii.PREFIX}/bii_2020.tif",
    "gcs_generation": "123",
    "gcs_crc32c": "AAAAAA==",
    "sha256": "f" * 64,
}
GRID = {
    "dtype": "float32",
    "nodata": None,
    "crs": "EPSG:4326",
    "transform": [0.0008983152841195215, 0.0, -180.0, 0.0, -0.0008983152841195215, 85.0],
    "shape": [161413, 400750],
}


def registered(manifest: dict) -> dict:
    return {
        "tilesets": [{"sources": [{"uris": [f"{OBJ['gcs_uri']}#123"]}]}],
        "properties": dict(manifest["properties"]),
        "startTime": manifest["startTime"],
        "endTime": manifest["endTime"],
        "bands": [
            {
                "id": "bii",
                "dataType": {"precision": "FLOAT"},
                "grid": {
                    "crsCode": "EPSG:4326",
                    "dimensions": {"width": 400750, "height": 161413},
                    "affineTransform": {
                        "scaleX": GRID["transform"][0],
                        "translateX": -180.0,
                        "scaleY": GRID["transform"][4],
                        "translateY": 85.0,
                    },
                },
            }
        ],
    }


def test_manifest_carries_credits_dates_and_unmodified_source_bands():
    manifest = bii.registration_manifest(2020, OBJ, GRID)
    assert manifest["name"] == "projects/landler-open-data/assets/datasets/bii/v1_1/bii_2020"
    assert (manifest["startTime"], manifest["endTime"]) == (
        "2020-01-01T00:00:00Z",
        "2021-01-01T00:00:00Z",
    )
    assert manifest["bands"] == [{"id": "bii", "tilesetId": "data", "tilesetBandIndex": 0}]
    properties = manifest["properties"]
    assert properties["license"] == "CC-BY-4.0"
    assert properties["citation"].startswith("Francis Gassert, Joe Mazzarello, and Samantha Hyde")
    assert properties["source_etag"] == "abc-2" and properties["gcs_generation"] == "123"
    assert all(k.replace("_", "").isalnum() for k in {**properties, **bii.collection_properties()})
    assert bii.CREDITS.items() <= bii.collection_properties().items()
    with_nodata = bii.registration_manifest(2020, OBJ, {**GRID, "nodata": 0.0})
    assert with_nodata["bands"][0]["missingData"] == {"values": [0.0]}


def test_registration_errors_detect_each_mismatch():
    manifest = bii.registration_manifest(2020, OBJ, GRID)
    assert bii.registration_errors(registered(manifest), manifest, OBJ, GRID) == []
    asset = registered(manifest)
    asset["tilesets"][0]["sources"][0]["uris"] = [OBJ["gcs_uri"] + "#999"]
    asset["properties"]["citation"] = "other"
    asset["bands"][0]["grid"]["affineTransform"]["translateY"] = 84.0
    asset["bands"][0]["grid"]["dimensions"]["width"] = 1
    errors = bii.registration_errors(asset, manifest, OBJ, GRID)
    assert len(errors) == 4


def test_choose_pixels_is_deterministic_and_covers_land_and_zero():
    overview = np.zeros((10, 20), dtype="float32")
    overview[2:5, 3:9] = 0.7
    pixels = bii.choose_pixels(overview, (1000, 2000), 2020, 20, 4)
    assert pixels == bii.choose_pixels(overview, (1000, 2000), 2020, 20, 4)
    assert len(pixels) == 20 and len(set(pixels)) == 20
    land = [p for p in pixels if overview[p[0] // 100, p[1] // 100] > 0]
    assert len(land) == 16
