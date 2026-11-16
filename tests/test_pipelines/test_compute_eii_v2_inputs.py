"""compute_eii_v2 reads annual NPP from the local CLMS DMP V2.0.2 sums."""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

rasterio = pytest.importorskip("rasterio")
from rasterio.transform import from_origin  # noqa: E402

sys.path.insert(0, str(Path(__file__).parents[2] / "pipelines" / "production"))
compute_eii_v2 = pytest.importorskip("compute_eii_v2")


@pytest.mark.parametrize("band", [1, 2])
def test_clms_vrt_exposes_one_band_of_the_annual_file(tmp_path, band):
    data = np.arange(64, dtype="int32").reshape(1, 8, 8)
    data = np.concatenate([data, data + 100])
    transform = from_origin(-180.0014880952381, 80.00148809523809, 1 / 336, 1 / 336)
    with rasterio.open(
        tmp_path / "annual_npp_2024.tif",
        "w",
        driver="GTiff",
        width=8,
        height=8,
        count=2,
        dtype="int32",
        crs="EPSG:4326",
        transform=transform,
        nodata=-1,
    ) as dst:
        dst.write(data)
    vrt = compute_eii_v2.clms_vrt(2024, tmp_path / "vrt" / "x.vrt", tmp_path, band=band)
    with rasterio.open(vrt) as src:
        assert src.count == 1
        assert src.transform == transform
        assert src.nodata == -1
        assert np.array_equal(src.read(1), data[band - 1])


def test_gcs_env_bills_the_user_project_and_uses_gdal_credentials(tmp_path, monkeypatch):
    gdal = ("GDAL_DISABLE_READDIR_ON_OPEN", "GDAL_HTTP_MAX_RETRY", "GDAL_HTTP_RETRY_DELAY")
    for key in ("GOOGLE_APPLICATION_CREDENTIALS", "GS_USER_PROJECT", "EII_BILLING_PROJECT", *gdal):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    with pytest.raises(RuntimeError, match="EII_BILLING_PROJECT"):
        compute_eii_v2.gcs_env()
    monkeypatch.setenv("EII_BILLING_PROJECT", "my-project")
    compute_eii_v2.gcs_env()
    assert os.environ["GS_USER_PROJECT"] == "my-project"
    assert "GOOGLE_APPLICATION_CREDENTIALS" not in os.environ
    adc = tmp_path / ".config/gcloud/application_default_credentials.json"
    adc.parent.mkdir(parents=True)
    adc.write_text('{"type": "service_account"}')
    compute_eii_v2.gcs_env()
    assert os.environ["GOOGLE_APPLICATION_CREDENTIALS"] == str(adc)
