"""Regenerate manuscript figures that are not emitted by analysis notebooks."""

from generate_abi_workflow import main as generate_abi_workflow
from generate_boundary_probability_residual_contrast import (
    main as generate_boundary_probability_residual_contrast,
)
from generate_deployment_support_training_distribution import (
    main as generate_deployment_support_training_distribution,
)


def main() -> None:
    generate_abi_workflow()
    generate_boundary_probability_residual_contrast()
    generate_deployment_support_training_distribution()


if __name__ == "__main__":
    main()
