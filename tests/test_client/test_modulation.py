"""Tests for eii.client.modulation module."""

import pytest

from eii.client.modulation import get_modulated_eii, get_nc_score
from eii.compute.modulation import apply_modulation, calculate_nc_score


class TestModulationClient:
    """Tests for modulation client functions."""

    def test_get_modulated_eii_scalar_kpis(self, amazon_small_polygon):
        """Validate scalar KPI modulation with a single GEE call."""
        kpis = {"biodiversity": 0.6, "soil": 0.4, "water": 0.8}

        result = get_modulated_eii(
            amazon_small_polygon,
            kpis=kpis,
            include_components=False,
            compute_mode="precomputed",
            scale=1000,
            stats=["mean"],
        )

        values = result["values"]
        base_mean = values["eii"]["mean"]
        assert isinstance(base_mean, int | float)

        expected_nc = calculate_nc_score(
            biodiversity_kpi=kpis["biodiversity"],
            soil_kpi=kpis["soil"],
            water_kpi=kpis["water"],
        )
        expected_modulated = apply_modulation(base_mean, expected_nc)

        assert values["nc_score"]["mean"] == pytest.approx(expected_nc)
        assert values["eii_modulated"]["mean"] == pytest.approx(expected_modulated)
        assert values["eii_modulated"]["mean"] >= base_mean
        assert 0 <= values["eii_modulated"]["mean"] <= 1

        assert values["biodiversity_kpi"]["mean"] == pytest.approx(kpis["biodiversity"])
        assert values["soil_kpi"]["mean"] == pytest.approx(kpis["soil"])
        assert values["water_kpi"]["mean"] == pytest.approx(kpis["water"])

    def test_modulation_input_validation(self, amazon_small_polygon):
        """Validate input errors without extra GEE calls."""
        with pytest.raises(ValueError, match="Must provide one of"):
            get_modulated_eii(amazon_small_polygon)

        with pytest.raises(ValueError, match="Provide only one of"):
            get_modulated_eii(
                amazon_small_polygon,
                kpis={"biodiversity": 0.5, "soil": 0.5, "water": 0.5},
                compute_default_kpis=True,
            )

        with pytest.raises(ValueError, match="Missing required KPI keys"):
            get_nc_score(amazon_small_polygon, kpis={"biodiversity": 0.5, "soil": 0.4})

        with pytest.raises(ValueError, match="must be numeric"):
            get_nc_score(
                amazon_small_polygon,
                kpis={"biodiversity": "bad", "soil": 0.4, "water": 0.5},
            )

        with pytest.raises(ValueError, match="must be in"):
            get_nc_score(
                amazon_small_polygon,
                kpis={"biodiversity": 1.2, "soil": 0.4, "water": 0.5},
            )

        with pytest.raises(ValueError, match="Weights must sum to 1.0"):
            get_nc_score(
                amazon_small_polygon,
                kpis={"biodiversity": 0.2, "soil": 0.3, "water": 0.4},
                kpi_weights={"biodiversity": 0.5, "soil": 0.5, "water": 0.5},
            )

        custom_weights = {"biodiversity": 0.2, "soil": 0.3, "water": 0.5}
        expected = calculate_nc_score(0.2, 0.3, 0.4, weights=custom_weights)
        assert get_nc_score(
            amazon_small_polygon,
            kpis={"biodiversity": 0.2, "soil": 0.3, "water": 0.4},
            kpi_weights=custom_weights,
        ) == pytest.approx(expected)
