"""Tests for eii.client.assets module."""

import pytest

from eii.client.assets import ASSETS, get_asset_info, get_asset_path


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

        structural_path = get_asset_path("structural_core_area")
        assert "structural" in structural_path.lower()

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
        required_assets = ["eii", "components", "structural_core_area"]

        for key in required_assets:
            assert key in ASSETS, f"Missing required asset: {key}"

        for name, info in ASSETS.items():
            assert "path" in info, f"Asset {name} missing 'path'"
            assert "description" in info, f"Asset {name} missing 'description'"
            # All defined assets should have retrievable paths
            assert isinstance(get_asset_path(name), str)
