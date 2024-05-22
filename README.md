# SCORT-ML-Pipeline
## Overview
This directory contains the code for a machine learning prediction model developed to identify colorectal cancer using a specific genetic signature. The signature was extracted using differential gene expression analysis tools such as Limma and samr. The model employs a gradient boosting machine trained on gene expression data from the S:CORT colorectal cancer cohort.

## Features
- **Prediction Model**: Utilizes gradient boosting to predict colorectal cancer based on gene expression data.
- **Signature Extraction**: Differential gene expression analysis tools like Limma and samr were used to identify the genetic signature.
- **R Shiny App**: The prediction tool is encapsulated in an R Shiny app, which can be run in RStudio or on an R Shiny server.

## Requirements
- **R**: Ensure you have R installed on your system.
- **RStudio**: Recommended for running the R Shiny app locally.
- **R Shiny Server**: Optionally, you can deploy the app on an R Shiny server for broader accessibility.
- **Gene Expression and Clinical Data**: Required for making predictions on validation sets.

## Usage
1. **Prepare Data**: Gather gene expression and clinical data for the study. Sample data can be found in the `Data/Validation` directory.
2. **Run the App**:
   - Open RStudio.
   - Load the R Shiny app script.
   - Run the app using the `runApp()` function.
3. **Make Predictions**: Use the app interface to input your validation set data and generate predictions.

## Directory Structure
- **`Data/Validation`**: Contains sample gene expression and clinical data for a validation set.

## Testing the Predictor
To test the predictor:
1. Start the R Shiny app.
2. Use the provided interface to load your validation set data.
3. Run predictions and view the results within the app.

## Contact
For any issues or questions, please contact [sr952@cantab.ac.uk].
