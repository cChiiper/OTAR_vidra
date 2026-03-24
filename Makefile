clean:
	rm -rf .nextflow* work results .pytest_cache
	find . -type d -name "__pycache__" -exec rm -rf {} +

docker_build:
	docker build --no-cache --platform=linux/amd64 -t headtext:latest .

docker_pull:
	docker pull ghcr.io/cchiiper/otar_vidra:latest