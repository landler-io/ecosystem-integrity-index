"""Read and write the model card of the potential-NPP model."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from importlib.resources import files
from pathlib import Path
from typing import Any


@dataclass
class ModelCard:
    """Provenance of the potential-NPP model.

    Predictors in inference order with their sources, hyperparameters, training data and
    cross-validated performance (``validation``).
    """

    model_id: str
    version: str
    response: str
    response_units: str
    predictors: list[str]
    predictor_sources: dict[str, str]
    hyperparameters: dict[str, Any]
    training_asset: str | None
    n_train: int | None
    sampling_design: dict[str, Any]
    validation: dict[str, Any]
    created: str | None
    notes: list[str]
    companion_models: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Reject incomplete identities and ambiguous predictor specifications."""
        for name in ("model_id", "version", "response", "response_units"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value:
                raise ValueError(f"{name} must be a nonempty string")
        if (
            not isinstance(self.predictors, list)
            or not self.predictors
            or any(not isinstance(name, str) or not name for name in self.predictors)
            or len(set(self.predictors)) != len(self.predictors)
        ):
            raise ValueError("predictors must be a nonempty ordered list of unique names")
        if not isinstance(self.predictor_sources, dict) or set(self.predictor_sources) != set(
            self.predictors
        ):
            raise ValueError("predictor_sources must cover exactly the model predictors")
        if any(
            not isinstance(source, str) or not source for source in self.predictor_sources.values()
        ):
            raise ValueError("predictor_sources must contain nonempty source descriptions")
        if self.n_train is not None and (type(self.n_train) is not int or self.n_train < 0):
            raise ValueError("n_train must be a nonnegative integer or null")

    def to_json(self) -> str:
        """Serialize the card to JSON."""
        return json.dumps(asdict(self), indent=2, allow_nan=False) + "\n"

    @classmethod
    def from_json(cls, value: str) -> ModelCard:
        """Parse a serialized model card."""
        return cls(**json.loads(value))

    def write(self, path: str | Path) -> None:
        """Write the card without creating or modifying model assets."""
        Path(path).write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def read(cls, path: str | Path) -> ModelCard:
        """Read a card from a user-provided JSON path."""
        return cls.from_json(Path(path).read_text(encoding="utf-8"))


def load_model_card(version: str = "v2") -> ModelCard:
    """Load a bundled card from a source checkout or installed wheel."""
    return ModelCard.from_json(
        files("eii.training")
        .joinpath("model_cards")
        .joinpath(f"{version}.json")
        .read_text(encoding="utf-8")
    )
