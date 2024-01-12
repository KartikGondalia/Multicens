# MultiCens

https://www.birdslab.in

## Getting Started

To use MultiCens, follow these steps:

1. **Download the GTEx Dataset**:
   - Visit [https://www.gtexportal.org/home/downloads](https://www.gtexportal.org/home/downloads) to download the GTEx dataset.
   - Save the downloaded dataset files into the `dataset/GTEx_Analysis_v8_eQTL_expression_matrices` folder within the project directory.

2. **Setting Up the Environment**:
   - Open the project directory in your terminal or command prompt.

   - Create a virtual environment for this project:
    
     python -m venv venv
 

   - Activate the virtual environment:
     - On Windows:
     
       venv\Scripts\activate
  

     - On macOS and Linux:
       
       source venv/bin/activate
  

3. **Installing Dependencies**:
   - Run the following command to install all the required libraries:
   
     pip install -r requirements.txt
     

   - Install the Ensembl version required for this project using the following command:
     
     pyensembl install --release 86 --species homo_sapiens
    

4. **Running the Code**:
   - Once you have set up the environment and installed the dependencies, you can run the MultiCens code using the following command:
  
     python app.py
     

## Backend Code
If you need to access the backend code, you can find it in the `util` folder. Inside the `util` folder, you will find two backend files:

- `util.py`: Contains utility functions and helper methods.
- `centrality.py`: Includes code related to centrality analysis.

## Disclaimer

This project and its code are the intellectual property of BirdsLab at IITM, and all rights are reserved. You may freely access and use the project for research and educational purposes. However, any commercial use or redistribution of this project without the explicit permission of BirdsLab is prohibited.


We hope MultiCens proves to be a valuable resource for your gene expression multicenter analysis. Thank you for using our tool!
