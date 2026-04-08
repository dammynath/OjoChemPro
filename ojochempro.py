"""
Ojo Chem Pro - Advanced 3D Chemistry Lab with Symmetry Analysis
Streamlit Deployment Version - Fixed for Python 3.13
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from stmol import showmol
import py3Dmol
import requests
import re
from datetime import datetime

# Page configuration
st.set_page_config(
    page_title="Ojo Chem Pro",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    /* Main container styling */
    .main-header {
        background: linear-gradient(135deg, #2c3e50 0%, #3498db 100%);
        padding: 1rem;
        border-radius: 10px;
        margin-bottom: 1rem;
    }
    
    .converter-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 1rem;
        border-radius: 10px;
        color: white;
        margin-bottom: 1rem;
    }
    
    .result-card {
        background: #e8f4f8;
        padding: 1rem;
        border-radius: 10px;
        margin-top: 1rem;
        border-left: 4px solid #667eea;
    }
    
    .molecule-card {
        background: white;
        padding: 0.75rem;
        border-radius: 10px;
        text-align: center;
        cursor: pointer;
        transition: all 0.2s;
        border: 1px solid #e0e0e0;
        margin-bottom: 0.5rem;
    }
    
    .molecule-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        border-color: #667eea;
    }
    
    .stButton button {
        width: 100%;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
    }
    
    .stButton button:hover {
        transform: scale(1.02);
    }
    
    /* Chat container styling */
    .chat-message-user {
        text-align: right;
        margin-bottom: 0.75rem;
    }
    
    .chat-message-user span {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 0.5rem 1rem;
        border-radius: 18px;
        display: inline-block;
        max-width: 90%;
        word-wrap: break-word;
    }
    
    .chat-message-bot {
        text-align: left;
        margin-bottom: 0.75rem;
    }
    
    .chat-message-bot span {
        background: #f0f2f6;
        color: #1e1e2f;
        padding: 0.5rem 1rem;
        border-radius: 18px;
        display: inline-block;
        max-width: 90%;
        border: 1px solid #e0e0e0;
        word-wrap: break-word;
    }
</style>
""", unsafe_allow_html=True)

# ========== COMPREHENSIVE PERIODIC TABLE DATABASE ==========
PERIODIC_TABLE = {
    "H": ["H", "Hydrogen", 1.008], "He": ["He", "Helium", 4.0026],
    "Li": ["Li", "Lithium", 6.94], "Be": ["Be", "Beryllium", 9.0122],
    "B": ["B", "Boron", 10.81], "C": ["C", "Carbon", 12.011],
    "N": ["N", "Nitrogen", 14.007], "O": ["O", "Oxygen", 15.999],
    "F": ["F", "Fluorine", 18.998], "Ne": ["Ne", "Neon", 20.18],
    "Na": ["Na", "Sodium", 22.99], "Mg": ["Mg", "Magnesium", 24.305],
    "Al": ["Al", "Aluminum", 26.982], "Si": ["Si", "Silicon", 28.086],
    "P": ["P", "Phosphorus", 30.974], "S": ["S", "Sulfur", 32.06],
    "Cl": ["Cl", "Chlorine", 35.45], "Ar": ["Ar", "Argon", 39.95],
    "K": ["K", "Potassium", 39.098], "Ca": ["Ca", "Calcium", 40.078],
    "Fe": ["Fe", "Iron", 55.845], "Cu": ["Cu", "Copper", 63.546],
    "Zn": ["Zn", "Zinc", 65.38], "Ag": ["Ag", "Silver", 107.87],
    "Au": ["Au", "Gold", 196.97], "Hg": ["Hg", "Mercury", 200.59],
    "Pb": ["Pb", "Lead", 207.2]
}

# Common compounds database
COMMON_COMPOUNDS = {
    "Water": {"formula": "H2O", "mass": 18.015},
    "Carbon Dioxide": {"formula": "CO2", "mass": 44.01},
    "Sodium Chloride": {"formula": "NaCl", "mass": 58.44},
    "Glucose": {"formula": "C6H12O6", "mass": 180.156},
    "Sulfuric Acid": {"formula": "H2SO4", "mass": 98.079},
    "Ammonia": {"formula": "NH3", "mass": 17.031},
    "Methane": {"formula": "CH4", "mass": 16.04},
    "Ethanol": {"formula": "C2H5OH", "mass": 46.07}
}

# Atomic masses dictionary
ATOMIC_MASSES = {symbol: data[2] for symbol, data in PERIODIC_TABLE.items()}

# Molecule symmetry database
MOLECULE_SYMMETRY = {
    "water": {"name": "Water", "formula": "H₂O", "point_group": "C₂v",
              "elements": [{"type": "C₂", "axis": [0, 1, 0], "description": "2-fold rotation axis"}],
              "description": "Bent molecule with C₂v symmetry."},
    "benzene": {"name": "Benzene", "formula": "C₆H₆", "point_group": "D₆h",
                "elements": [{"type": "C₆", "axis": [0, 1, 0], "description": "6-fold principal axis"}],
                "description": "Aromatic ring with D₆h symmetry."},
    "methane": {"name": "Methane", "formula": "CH₄", "point_group": "Td",
                "elements": [{"type": "C₃", "axis": [1, 1, 1], "description": "3-fold axes"}],
                "description": "Perfect tetrahedral molecule with Td symmetry."}
}

# Molecule 3D structures
MOLECULE_STRUCTURES = {
    "water": "O", "benzene": "c1ccccc1", "methane": "C",
    "ammonia": "N", "carbon_dioxide": "O=C=O", "ethanol": "CCO"
}

# ========== UTILITY FUNCTIONS ==========
def format_formula(formula):
    """Convert chemical formula to subscript format"""
    subscripts = {'0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄',
                  '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉'}
    for num, sub in subscripts.items():
        formula = formula.replace(num, sub)
    return formula

def calculate_molar_mass(formula):
    """Calculate molar mass from chemical formula"""
    try:
        total_mass = 0
        i = 0
        elements = []
        
        while i < len(formula):
            if i + 1 < len(formula) and formula[i + 1].islower():
                element = formula[i:i + 2]
                i += 2
            else:
                element = formula[i]
                i += 1
            
            count = ''
            while i < len(formula) and formula[i].isdigit():
                count += formula[i]
                i += 1
            
            num = int(count) if count else 1
            
            if element in ATOMIC_MASSES:
                mass = ATOMIC_MASSES[element] * num
                total_mass += mass
                elements.append(f"{element}{num if num > 1 else ''}: {mass:.2f} g/mol")
            else:
                return None, f"Unknown element: {element}"
        
        return total_mass, elements
    except Exception as e:
        return None, str(e)

def mass_to_moles(mass, molar_mass):
    """Convert mass to moles"""
    if molar_mass > 0:
        return mass / molar_mass
    return None

def moles_to_mass(moles, molar_mass):
    """Convert moles to mass"""
    return moles * molar_mass

def balance_equation(equation):
    """Balance chemical equations"""
    balanced_equations = {
        "h2+o2->h2o": "2H₂ + O₂ → 2H₂O",
        "ch4+o2->co2+h2o": "CH₄ + 2O₂ → CO₂ + 2H₂O",
        "na+cl2->nacl": "2Na + Cl₂ → 2NaCl",
        "n2+h2->nh3": "N₂ + 3H₂ → 2NH₃"
    }
    key = equation.lower().replace(" ", "").replace("->", "->")
    return balanced_equations.get(key, None)

def search_pubchem(compound):
    """Search PubChem API"""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound}/property/MolecularFormula,MolecularWeight/JSON"
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
                prop = data["PropertyTable"]["Properties"][0]
                return {"formula": prop.get("MolecularFormula", "N/A"),
                        "weight": prop.get("MolecularWeight", "N/A")}
    except:
        pass
    return None

def get_chatbot_response(query, safe_mode=True):
    """Generate response for user queries"""
    query_lower = query.lower()
    
    hazardous_keywords = ['explosive', 'cyanide', 'mustard gas', 'weapon', 'toxic gas']
    if safe_mode and any(keyword in query_lower for keyword in hazardous_keywords):
        return "⚠️ **Safety Notice**: This topic involves hazardous materials and is blocked in Safe Mode."
    
    for molecule in MOLECULE_SYMMETRY.keys():
        if molecule in query_lower:
            sym = MOLECULE_SYMMETRY[molecule]
            return f"🔬 **{sym['name']}**\n\n**Formula:** {sym['formula']}\n**Point Group:** {sym['point_group']}\n\n{sym['description']}"
    
    if "balance" in query_lower or "->" in query:
        balanced = balance_equation(query)
        if balanced:
            return f"✅ **Balanced Equation:**\n\n{balanced}"
        return "🤔 Try: H2 + O2 -> H2O or CH4 + O2 -> CO2 + H2O"
    
    if "molar mass" in query_lower:
        formula_match = re.search(r'([A-Z][a-z]?\d*)+', query.upper())
        if formula_match:
            formula = formula_match.group()
            mass, elements = calculate_molar_mass(formula)
            if mass:
                return f"📊 **Molar Mass of {format_formula(formula)}**\n\n**Molar Mass:** {mass:.4f} g/mol"
        return "Please provide a chemical formula (e.g., 'Molar mass of H2O')"
    
    if "pubchem" in query_lower:
        compound = query_lower.replace("pubchem", "").strip()
        if compound:
            result = search_pubchem(compound)
            if result:
                return f"🔍 **PubChem Results for {compound}**\n\n📝 Formula: {result['formula']}\n⚖️ Weight: {result['weight']} g/mol"
        return "Please specify a compound to search"
    
    return f"🧪 **Chemistry Assistant**\n\nI can help with:\n• 3D Molecule Viewer\n• Symmetry Analysis\n• Mass-Mole Conversions\n• Equation Balancing\n• Molar Mass Calculator\n\n**Try:** 'Show me benzene' or 'Balance: H2 + O2 -> H2O'"

# ========== 3D MOLECULE VISUALIZATION ==========
def create_3d_molecule(smiles, molecule_name):
    """Create 3D molecule visualization using py3Dmol"""
    viewer = py3Dmol.view(width=500, height=400)
    viewer.addModel(smiles, 'smi')
    viewer.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    viewer.setBackgroundColor('0x1a1a2e')
    viewer.zoomTo()
    viewer.spin(True)
    return viewer

def create_symmetry_visualization(molecule_name):
    """Create symmetry element visualization"""
    sym_data = MOLECULE_SYMMETRY.get(molecule_name.lower())
    if not sym_data:
        return None
    
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(x=[0], y=[0], z=[0], mode='markers',
                                marker=dict(size=8, color='white'), name='Center'))
    
    for element in sym_data["elements"]:
        if "axis" in element:
            axis = element["axis"]
            fig.add_trace(go.Scatter3d(x=[-axis[0]*2, axis[0]*2],
                                        y=[-axis[1]*2, axis[1]*2],
                                        z=[-axis[2]*2, axis[2]*2],
                                        mode='lines', line=dict(color='red', width=4),
                                        name=f"{element['type']} Axis"))
    
    fig.update_layout(title=f"Symmetry of {sym_data['name']} ({sym_data['point_group']})",
                      scene=dict(bgcolor='#1a1a2e', xaxis=dict(color='white'),
                                 yaxis=dict(color='white'), zaxis=dict(color='white')),
                      paper_bgcolor='#1a1a2e', font=dict(color='white'), height=450)
    return fig

# ========== MAIN APPLICATION ==========
def main():
    # Header
    st.markdown("""
    <div class="main-header">
        <h1 style="color: white; margin: 0;">🧪 Ojo Chem Pro</h1>
        <p style="color: white; margin: 0; opacity: 0.9;">Advanced 3D Chemistry Lab with Mass-Mole Converter</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Emergency banner (minimized)
    with st.expander("🚨 Emergency Information", expanded=False):
        st.warning("**Poison Control:** 1-800-222-1222 | **CHEMTREC:** 1-800-424-9300\n\n**First Aid:** Flush with water for 15 minutes, seek immediate medical attention.")
    
    # Sidebar for mode and mass-mole converter
    with st.sidebar:
        st.markdown("### 🎮 Controls")
        safe_mode = st.toggle("🔒 Safe Mode", value=True)
        if not safe_mode:
            st.warning("⚠️ Expert Mode: Adult supervision recommended.")
        
        st.divider()
        
        # Mass to Mole Converter
        st.markdown("### ⚖️ Mass ↔ Mole Converter")
        
        conversion_type = st.radio("Select conversion:", ["Mass → Moles", "Moles → Mass"], horizontal=True)
        
        # Compound selection
        compound_option = st.radio("Compound:", ["Common", "Custom Formula"], horizontal=True)
        
        molar_mass = None
        formula_display = ""
        
        if compound_option == "Common":
            compound_name = st.selectbox("Select:", list(COMMON_COMPOUNDS.keys()))
            formula_display = COMMON_COMPOUNDS[compound_name]["formula"]
            molar_mass = COMMON_COMPOUNDS[compound_name]["mass"]
            st.info(f"**Formula:** {format_formula(formula_display)} | **Molar Mass:** {molar_mass:.3f} g/mol")
        else:
            formula_input = st.text_input("Formula:", placeholder="e.g., H2O, C6H12O6")
            if formula_input:
                mass_calc, _ = calculate_molar_mass(formula_input.upper())
                if mass_calc:
                    molar_mass = mass_calc
                    formula_display = formula_input
                    st.success(f"**Molar Mass:** {molar_mass:.3f} g/mol")
        
        # Conversion calculator
        if molar_mass and molar_mass > 0:
            st.divider()
            if conversion_type == "Mass → Moles":
                mass_input = st.number_input("Mass (grams):", min_value=0.0, value=1.0, step=0.1)
                if mass_input > 0:
                    moles = mass_to_moles(mass_input, molar_mass)
                    st.markdown(f"""
                    <div class="result-card">
                        <strong>📊 Result:</strong><br>
                        Mass: {mass_input:.4f} g<br>
                        Molar Mass: {molar_mass:.4f} g/mol<br>
                        <strong>Moles = {moles:.6f} mol</strong><br>
                        Molecules: {moles * 6.022e23:.2e}
                    </div>
                    """, unsafe_allow_html=True)
            else:
                moles_input = st.number_input("Moles:", min_value=0.0, value=1.0, step=0.1)
                if moles_input > 0:
                    mass = moles_to_mass(moles_input, molar_mass)
                    st.markdown(f"""
                    <div class="result-card">
                        <strong>📊 Result:</strong><br>
                        Moles: {moles_input:.4f} mol<br>
                        Molar Mass: {molar_mass:.4f} g/mol<br>
                        <strong>Mass = {mass:.4f} g</strong>
                    </div>
                    """, unsafe_allow_html=True)
        
        st.divider()
        
        # Quick molecules
        st.markdown("### 📚 Quick Molecules")
        molecules = ["Water", "Benzene", "Methane", "Ammonia", "CO₂", "Ethanol"]
        cols = st.columns(3)
        for i, mol in enumerate(molecules):
            with cols[i % 3]:
                if st.button(f"🔬 {mol}", key=f"quick_{mol}"):
                    st.session_state['selected_molecule'] = mol.lower()
                    st.rerun()
        
        st.divider()
        st.markdown(f"**Elements:** {len(PERIODIC_TABLE)} | **Compounds:** {len(COMMON_COMPOUNDS)}")
    
    # Main content area - Two columns for content and chat
    col_left, col_right = st.columns([2.2, 1.2])
    
    # LEFT COLUMN - Main content
    with col_left:
        tab1, tab2, tab3 = st.tabs(["🔬 3D Viewer", "🔄 Symmetry Analysis", "⚙️ Tools"])
        
        with tab1:
            molecule_name = st.selectbox("Select Molecule:", list(MOLECULE_STRUCTURES.keys()),
                                          format_func=lambda x: x.title())
            if 'selected_molecule' in st.session_state:
                molecule_name = st.session_state['selected_molecule']
            
            smiles = MOLECULE_STRUCTURES.get(molecule_name, "C")
            viewer = create_3d_molecule(smiles, molecule_name)
            showmol(viewer, height=450, width=650)
            
            if molecule_name in MOLECULE_SYMMETRY:
                sym = MOLECULE_SYMMETRY[molecule_name]
                st.info(f"**{sym['name']}** ({sym['formula']}) - **Point Group:** {sym['point_group']}")
        
        with tab2:
            sym_molecule = st.selectbox("Select for Analysis:", list(MOLECULE_SYMMETRY.keys()),
                                         format_func=lambda x: MOLECULE_SYMMETRY[x]['name'])
            fig = create_symmetry_visualization(sym_molecule)
            if fig:
                st.plotly_chart(fig, use_container_width=True)
            
            sym_data = MOLECULE_SYMMETRY.get(sym_molecule, {})
            st.markdown(f"**Point Group:** `{sym_data.get('point_group', 'N/A')}`")
            for elem in sym_data.get('elements', []):
                st.markdown(f"- **{elem['type']}**: {elem['description']}")
        
        with tab3:
            tool_tabs = st.tabs(["📊 Molar Mass", "⚖️ Equation Balancer", "🔍 PubChem"])
            
            with tool_tabs[0]:
                formula_input = st.text_input("Formula:", placeholder="H2O, C6H12O6, NaCl")
                if st.button("Calculate Molar Mass"):
                    if formula_input:
                        mass, elements = calculate_molar_mass(formula_input.upper())
                        if mass:
                            st.success(f"**Molar Mass:** {mass:.4f} g/mol")
                            with st.expander("Composition"):
                                for elem in elements:
                                    st.write(elem)
            
            with tool_tabs[1]:
                eq_input = st.text_input("Equation:", placeholder="H2 + O2 -> H2O")
                if st.button("Balance"):
                    balanced = balance_equation(eq_input)
                    if balanced:
                        st.success(f"**Balanced:** {balanced}")
                    else:
                        st.info("Try: H2 + O2 -> H2O or CH4 + O2 -> CO2 + H2O")
            
            with tool_tabs[2]:
                compound_search = st.text_input("Compound:", placeholder="aspirin, caffeine, glucose")
                if st.button("Search PubChem"):
                    if compound_search:
                        with st.spinner("Searching..."):
                            result = search_pubchem(compound_search)
                            if result:
                                st.success(f"**Formula:** {result['formula']}\n**Weight:** {result['weight']} g/mol")
    
    # RIGHT COLUMN - Chat interface (NO st.chat_input inside column)
    with col_right:
        st.markdown("### 💬 Chemistry Assistant")
        
        # Initialize chat history
        if "chat_messages" not in st.session_state:
            st.session_state.chat_messages = [
                {"role": "assistant", "content": "👋 Hi! Ask me about molecules, symmetry, or chemistry!"}
            ]
        
        # Display chat messages
        chat_container = st.container()
        with chat_container:
            st.markdown('<div style="height: 450px; overflow-y: auto; padding: 1rem; background: #f8f9fa; border-radius: 10px; margin-bottom: 1rem;">', unsafe_allow_html=True)
            for msg in st.session_state.chat_messages:
                if msg["role"] == "user":
                    st.markdown(f'<div class="chat-message-user"><span>{msg["content"]}</span></div>', unsafe_allow_html=True)
                else:
                    st.markdown(f'<div class="chat-message-bot"><span>{msg["content"]}</span></div>', unsafe_allow_html=True)
            st.markdown('</div>', unsafe_allow_html=True)
        
        # Quick suggestions
        st.markdown("#### Quick Suggestions")
        c1, c2 = st.columns(2)
        with c1:
            if st.button("🔬 Show benzene", use_container_width=True):
                st.session_state.chat_messages.append({"role": "user", "content": "Show me benzene"})
                response = get_chatbot_response("Show me benzene", safe_mode)
                st.session_state.chat_messages.append({"role": "assistant", "content": response})
                st.rerun()
            if st.button("⚖️ Balance H2+O2", use_container_width=True):
                st.session_state.chat_messages.append({"role": "user", "content": "Balance: H2 + O2 -> H2O"})
                response = get_chatbot_response("Balance: H2 + O2 -> H2O", safe_mode)
                st.session_state.chat_messages.append({"role": "assistant", "content": response})
                st.rerun()
        
        with c2:
            if st.button("🔄 Water symmetry", use_container_width=True):
                st.session_state.chat_messages.append({"role": "user", "content": "Symmetry of water"})
                response = get_chatbot_response("Symmetry of water", safe_mode)
                st.session_state.chat_messages.append({"role": "assistant", "content": response})
                st.rerun()
            if st.button("📊 Molar mass H2O", use_container_width=True):
                st.session_state.chat_messages.append({"role": "user", "content": "Molar mass of H2O"})
                response = get_chatbot_response("Molar mass of H2O", safe_mode)
                st.session_state.chat_messages.append({"role": "assistant", "content": response})
                st.rerun()
    
    # CHAT INPUT - Placed OUTSIDE all columns, at the bottom of the app
    st.divider()
    
    # Create a container for chat input at the bottom
    chat_input_container = st.container()
    with chat_input_container:
        user_input = st.chat_input("Ask me about chemistry...", key="chemistry_chat_input_bottom")
        
        if user_input:
            st.session_state.chat_messages.append({"role": "user", "content": user_input})
            with st.spinner("Thinking..."):
                response = get_chatbot_response(user_input, safe_mode)
            st.session_state.chat_messages.append({"role": "assistant", "content": response})
            st.rerun()

if __name__ == "__main__":
    main()
