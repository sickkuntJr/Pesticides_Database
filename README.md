# Pesticide Database Project

This project scrapes pesticide data from a website, generates molecular structure images, and displays the data in a GUI application.

The data is provided, but the scraping and image generation scripts are provided in src if you wish to do it yourself.

## Project Structure

pesticides_project/
│
├── data/
│ ├── pesticides_data.json # JSON data file
│ ├── molecule_images/ # Directory for generated images
│
├── src/
│ ├── scrape_pesticides.py # Script for web scraping
│ └── generate_molecules.py # Script for molecule generation
│
├─── app.py # Main application script
│
├── requirements.txt # List of dependencies
└── README.md # Project documentation



## Setup

1. **Create and activate a virtual environment**:

    ```bash
    conda create -n .venv
    conda activate .venv
    ```

2. **Install dependencies**:

    ```bash
    pip install -r requirements.txt
    ```

3. **Run the scraping script**:

    ```bash
    python src/scrape_pesticides.py
    ```

4. **Generate molecular structure images**:

    ```bash
    python src/generate_molecules.py
    ```

5. **Run the application**:

    ```bash
    python src/app.py
    ```
