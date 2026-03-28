# Create virtual environment
python -m venv chem_env
source chem_env/bin/activate  # On Windows: chem_env\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run the app
streamlit run chem_assistant.py