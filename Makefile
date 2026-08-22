.PHONY: setup lint test check build clean publish publish-test

setup:
	uv sync

lint:
	uvx ruff check

test:
	uv run pytest

check: lint test

build:
	uv build

clean:
	rm -rf dist

# Publish dist/* to PyPI. Requires an API token in the environment:
#   export UV_PUBLISH_TOKEN=pypi-<token>
publish: clean check build
	uv publish

# Publish dist/* to TestPyPI (use a token from test.pypi.org).
publish-test: clean check build
	uv publish --publish-url https://test.pypi.org/legacy/