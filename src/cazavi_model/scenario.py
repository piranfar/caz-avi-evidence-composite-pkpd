"""Validated metadata schema for reproducible simulation scenarios."""

from __future__ import annotations

from pydantic import BaseModel, Field, PositiveInt, field_validator


class TargetDefinition(BaseModel):
    component: str
    metric: str
    threshold: float
    unit: str | None = None


class ScenarioMetadata(BaseModel):
    scenario_id: str = Field(min_length=1)
    source_studies: list[str]
    virtual_population_size: PositiveInt
    random_seed: int = Field(ge=0)
    population_specification_version: str
    pk_model_version: str
    regimen_version: str
    mic_registry_version: str
    free_fraction_policy: str
    targets: list[TargetDefinition]
    toxicity_rule: str | None = None
    synthetic_assumptions: list[str] = []
    validation_benchmarks: list[str] = []

    @field_validator("source_studies")
    @classmethod
    def require_source_studies(cls, value: list[str]) -> list[str]:
        if not value:
            raise ValueError("At least one source study must be recorded.")
        return value
