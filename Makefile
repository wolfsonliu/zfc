.PHONY: setup lint test check

setup:
	uv sync

lint:
	uvx ruff check

test:
	uv run pytest

check: lint test
