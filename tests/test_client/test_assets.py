"""Tests for eii.client.assets module."""

import pytest

from eii.client.assets import (
    ASSETS,
    get_asset_info,
    get_asset_path,
    get_version,
    set_version,
    use_version,
)


class TestAssets:
    """Tests for asset path and info retrieval."""

    def test_asset_path_retrieval(self):
        """Test asset path retrieval for valid and invalid keys."""
        # Valid assets return string paths
        eii_path = get_asset_path("eii")
        assert isinstance(eii_path, str)
        assert eii_path.startswith("projects/")

        components_path = get_asset_path("components")
        assert isinstance(components_path, str)

        assert get_asset_path("eii", version="v2").endswith("products/eii/v2_0/release/global")
        assert "v2" not in get_asset_info("eii", version="v2")

        assert get_asset_path("natural_npp", version="v2").endswith("potential_npp_olci_gapfill")

        # Unknown asset raises KeyError with helpful message
        with pytest.raises(KeyError) as exc_info:
            get_asset_path("nonexistent_asset")
        assert "nonexistent_asset" in str(exc_info.value)
        assert "Available" in str(exc_info.value)

    def test_asset_info_retrieval(self):
        """Test asset info retrieval returns expected structure."""
        # EII info has required fields
        eii_info = get_asset_info("eii")
        assert "path" in eii_info
        assert "description" in eii_info
        assert "bands" in eii_info
        assert "eii" in eii_info["bands"]

        # Components info has component bands
        comp_info = get_asset_info("components")
        assert all(
            band in comp_info["bands"]
            for band in ["functional_integrity", "structural_integrity", "compositional_integrity"]
        )

        # Returns a copy, not original
        eii_info["test_key"] = "test_value"
        assert "test_key" not in get_asset_info("eii")

        # Unknown raises KeyError
        with pytest.raises(KeyError):
            get_asset_info("nonexistent")

    def test_assets_dict_completeness(self):
        """Test ASSETS dict has required structure."""
        required_assets = ["eii", "components", "wetland_flag", "natural_npp"]

        for key in required_assets:
            assert key in ASSETS, f"Missing required asset: {key}"

        for name, info in ASSETS.items():
            assert "path" in info, f"Asset {name} missing 'path'"
            assert "description" in info, f"Asset {name} missing 'description'"
            # All defined assets should have retrievable paths
            assert isinstance(get_asset_path(name), str)


class TestVersionSelection:
    """Tests for dataset version resolution."""

    @pytest.fixture(autouse=True)
    def _restore_version(self):
        """Ensure the session version is restored after each test."""
        original = get_version()
        yield
        set_version(original)

    def test_explicit_version_argument(self):
        """An explicit version is embedded in the path and beats the session default."""
        assert get_asset_path("eii", version="v1").endswith("eii_global_v1")
        assert get_asset_path("eii", version="v9").endswith("eii_global_v9")
        assert get_asset_path("natural_npp", version="v1").startswith(
            "projects/landler-open-data/assets/eii/intermediate/"
        )

    def test_explicit_version_does_not_change_session(self):
        """Passing version per call must not mutate the active session version."""
        before = get_version()
        get_asset_path("eii", version="v1")
        assert get_version() == before

    def test_default_version_is_v2(self):
        """Without config overrides or set_version, retrieval reads v2."""
        from eii.config import VERSION

        assert VERSION == "v2"
        assert get_asset_path("eii", version=VERSION).endswith("products/eii/v2_0/release/global")

    def test_set_version_changes_default(self):
        """set_version changes the version used when none is passed."""
        original = get_version()
        try:
            set_version("v1")
            assert get_version() == "v1"
            assert get_asset_path("eii").endswith("eii_global_v1")
        finally:
            set_version(original)

    def test_use_version_context_manager(self):
        """use_version applies within the block and restores afterwards."""
        original = get_version()
        with use_version("v1"):
            assert get_version() == "v1"
            assert get_asset_path("components").endswith("eii_global_v1")
        assert get_version() == original

    def test_use_version_restores_on_exception(self):
        """The session version is restored even if the block raises."""
        original = get_version()
        with pytest.raises(RuntimeError), use_version("v1"):
            raise RuntimeError("boom")
        assert get_version() == original

    def test_set_version_rejects_empty(self):
        """An empty or non-string version is rejected."""
        with pytest.raises(ValueError):
            set_version("")
