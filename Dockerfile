FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && \
    apt-get install -y \
    clustalo \
    procps \
    curl \
    ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Create directory structure
RUN mkdir -p /app/kinase_project/{data,results}

# Set default email (override with -e flag)
ENV ENTREZ_EMAIL="tahagill99@gmail.com"

# Start script with immediate timeout protection
CMD ["timeout", "1h", "python", "main.py"]