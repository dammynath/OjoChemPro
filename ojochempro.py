"""
Ojo Chem Pro - Advanced 3D Chemistry Lab with Mass-Mole Converter
Streamlit Deployment Version - Adaptive Theming & Complete Periodic Table
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

# ========== ADAPTIVE THEME DETECTION ==========
def get_theme_colors():
    """Get theme-adaptive colors based on Streamlit theme"""
    # Detect if dark mode is active
    is_dark = st.get_option("theme.base") == "dark"
    
    if is_dark:
        return {
            "bg_color": "#1e1e2f",
            "card_bg": "#2d2d3f",
            "text_color": "#ffffff",
            "text_secondary": "#b0b0b0",
            "border_color": "#404060",
            "result_bg": "#2d2d3f",
            "result_border": "#667eea",
            "success_bg": "#1a3a2a",
            "success_text": "#00ff88",
            "warning_bg": "#3a2a1a",
            "warning_text": "#ffaa44",
            "info_bg": "#1a2a3a",
            "info_text": "#44aaff"
        }
    else:
        return {
            "bg_color": "#ffffff",
            "card_bg": "#f8f9fa",
            "text_color": "#1e1e2f",
            "text_secondary": "#666666",
            "border_color": "#e0e0e0",
            "result_bg": "#e8f4f8",
            "result_border": "#667eea",
            "success_bg": "#d4edda",
            "success_text": "#155724",
            "warning_bg": "#fff3cd",
            "warning_text": "#856404",
            "info_bg": "#d1ecf1",
            "info_text": "#0c5460"
        }

# Custom CSS with theme variables
st.markdown("""
<style>
    /* Main container styling */
    .main-header {
        background: linear-gradient(135deg, #2c3e50 0%, #3498db 100%);
        padding: 1rem;
        border-radius: 10px;
        margin-bottom: 1rem;
    }
    
    /* Adaptive result card */
    .result-card {
        background: var(--result-bg, #e8f4f8);
        padding: 1rem;
        border-radius: 10px;
        margin-top: 1rem;
        border-left: 4px solid var(--result-border, #667eea);
        color: var(--text-color, #1e1e2f);
        transition: all 0.3s ease;
    }
    
    .result-card strong {
        color: var(--result-border, #667eea);
    }
    
    /* Adaptive converter card */
    .converter-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 1rem;
        border-radius: 10px;
        color: white;
        margin-bottom: 1rem;
    }
    
    /* Adaptive molecule card */
    .molecule-card {
        background: var(--card-bg, white);
        padding: 0.75rem;
        border-radius: 10px;
        text-align: center;
        cursor: pointer;
        transition: all 0.2s;
        border: 1px solid var(--border-color, #e0e0e0);
        margin-bottom: 0.5rem;
        color: var(--text-color, #1e1e2f);
    }
    
    .molecule-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        border-color: #667eea;
    }
    
    /* Chat message styling */
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
        background: var(--card-bg, #f0f2f6);
        color: var(--text-color, #1e1e2f);
        padding: 0.5rem 1rem;
        border-radius: 18px;
        display: inline-block;
        max-width: 90%;
        border: 1px solid var(--border-color, #e0e0e0);
        word-wrap: break-word;
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
    
    /* Periodic table styling */
    .periodic-grid {
        display: grid;
        grid-template-columns: repeat(18, 1fr);
        gap: 5px;
        margin: 10px 0;
    }
    
    .element-cell {
        background: var(--card-bg, #f8f9fa);
        border: 1px solid var(--border-color, #e0e0e0);
        border-radius: 5px;
        padding: 5px;
        text-align: center;
        cursor: pointer;
        transition: all 0.2s;
        font-size: 11px;
        color: var(--text-color, #1e1e2f);
    }
    
    .element-cell:hover {
        transform: scale(1.05);
        border-color: #667eea;
        box-shadow: 0 2px 8px rgba(102, 126, 234, 0.3);
    }
    
    .element-symbol {
        font-weight: bold;
        font-size: 14px;
    }
    
    .element-number {
        font-size: 9px;
        opacity: 0.7;
    }
    
    /* Lanthanide/Actinide series */
    .series-grid {
        display: grid;
        grid-template-columns: repeat(15, 1fr);
        gap: 5px;
        margin: 10px 0;
    }
</style>
""", unsafe_allow_html=True)

# ========== COMPLETE PERIODIC TABLE WITH ALL TRANSITION METALS ==========
PERIODIC_TABLE = {
    # Period 1
    "H": ["H", "Hydrogen", 1.008, 1, 1], "He": ["He", "Helium", 4.0026, 18, 1],
    # Period 2
    "Li": ["Li", "Lithium", 6.94, 1, 2], "Be": ["Be", "Beryllium", 9.0122, 2, 2],
    "B": ["B", "Boron", 10.81, 13, 2], "C": ["C", "Carbon", 12.011, 14, 2],
    "N": ["N", "Nitrogen", 14.007, 15, 2], "O": ["O", "Oxygen", 15.999, 16, 2],
    "F": ["F", "Fluorine", 18.998, 17, 2], "Ne": ["Ne", "Neon", 20.18, 18, 2],
    # Period 3
    "Na": ["Na", "Sodium", 22.99, 1, 3], "Mg": ["Mg", "Magnesium", 24.305, 2, 3],
    "Al": ["Al", "Aluminum", 26.982, 13, 3], "Si": ["Si", "Silicon", 28.086, 14, 3],
    "P": ["P", "Phosphorus", 30.974, 15, 3], "S": ["S", "Sulfur", 32.06, 16, 3],
    "Cl": ["Cl", "Chlorine", 35.45, 17, 3], "Ar": ["Ar", "Argon", 39.95, 18, 3],
    # Period 4 - Including Transition Metals
    "K": ["K", "Potassium", 39.098, 1, 4], "Ca": ["Ca", "Calcium", 40.078, 2, 4],
    "Sc": ["Sc", "Scandium", 44.956, 3, 4], "Ti": ["Ti", "Titanium", 47.867, 4, 4],
    "V": ["V", "Vanadium", 50.942, 5, 4], "Cr": ["Cr", "Chromium", 51.996, 6, 4],
    "Mn": ["Mn", "Manganese", 54.938, 7, 4], "Fe": ["Fe", "Iron", 55.845, 8, 4],
    "Co": ["Co", "Cobalt", 58.933, 9, 4], "Ni": ["Ni", "Nickel", 58.693, 10, 4],
    "Cu": ["Cu", "Copper", 63.546, 11, 4], "Zn": ["Zn", "Zinc", 65.38, 12, 4],
    "Ga": ["Ga", "Gallium", 69.723, 13, 4], "Ge": ["Ge", "Germanium", 72.63, 14, 4],
    "As": ["As", "Arsenic", 74.922, 15, 4], "Se": ["Se", "Selenium", 78.96, 16, 4],
    "Br": ["Br", "Bromine", 79.904, 17, 4], "Kr": ["Kr", "Krypton", 83.798, 18, 4],
    # Period 5 - Including Transition Metals
    "Rb": ["Rb", "Rubidium", 85.468, 1, 5], "Sr": ["Sr", "Strontium", 87.62, 2, 5],
    "Y": ["Y", "Yttrium", 88.906, 3, 5], "Zr": ["Zr", "Zirconium", 91.224, 4, 5],
    "Nb": ["Nb", "Niobium", 92.906, 5, 5], "Mo": ["Mo", "Molybdenum", 95.95, 6, 5],
    "Tc": ["Tc", "Technetium", 98, 7, 5], "Ru": ["Ru", "Ruthenium", 101.07, 8, 5],
    "Rh": ["Rh", "Rhodium", 102.91, 9, 5], "Pd": ["Pd", "Palladium", 106.42, 10, 5],
    "Ag": ["Ag", "Silver", 107.87, 11, 5], "Cd": ["Cd", "Cadmium", 112.41, 12, 5],
    "In": ["In", "Indium", 114.82, 13, 5], "Sn": ["Sn", "Tin", 118.71, 14, 5],
    "Sb": ["Sb", "Antimony", 121.76, 15, 5], "Te": ["Te", "Tellurium", 127.6, 16, 5],
    "I": ["I", "Iodine", 126.90, 17, 5], "Xe": ["Xe", "Xenon", 131.29, 18, 5],
    # Period 6 - Including Lanthanides and Transition Metals
    "Cs": ["Cs", "Cesium", 132.91, 1, 6], "Ba": ["Ba", "Barium", 137.33, 2, 6],
    "La": ["La", "Lanthanum", 138.91, 3, 6], "Ce": ["Ce", "Cerium", 140.12, 101, 6],
    "Pr": ["Pr", "Praseodymium", 140.91, 101, 6], "Nd": ["Nd", "Neodymium", 144.24, 101, 6],
    "Pm": ["Pm", "Promethium", 145, 101, 6], "Sm": ["Sm", "Samarium", 150.36, 101, 6],
    "Eu": ["Eu", "Europium", 151.96, 101, 6], "Gd": ["Gd", "Gadolinium", 157.25, 101, 6],
    "Tb": ["Tb", "Terbium", 158.93, 101, 6], "Dy": ["Dy", "Dysprosium", 162.5, 101, 6],
    "Ho": ["Ho", "Holmium", 164.93, 101, 6], "Er": ["Er", "Erbium", 167.26, 101, 6],
    "Tm": ["Tm", "Thulium", 168.93, 101, 6], "Yb": ["Yb", "Ytterbium", 173.05, 101, 6],
    "Lu": ["Lu", "Lutetium", 174.97, 3, 6],
    "Hf": ["Hf", "Hafnium", 178.49, 4, 6], "Ta": ["Ta", "Tantalum", 180.95, 5, 6],
    "W": ["W", "Tungsten", 183.84, 6, 6], "Re": ["Re", "Rhenium", 186.21, 7, 6],
    "Os": ["Os", "Osmium", 190.23, 8, 6], "Ir": ["Ir", "Iridium", 192.22, 9, 6],
    "Pt": ["Pt", "Platinum", 195.08, 10, 6], "Au": ["Au", "Gold", 196.97, 11, 6],
    "Hg": ["Hg", "Mercury", 200.59, 12, 6], "Tl": ["Tl", "Thallium", 204.38, 13, 6],
    "Pb": ["Pb", "Lead", 207.2, 14, 6], "Bi": ["Bi", "Bismuth", 208.98, 15, 6],
    "Po": ["Po", "Polonium", 209, 16, 6], "At": ["At", "Astatine", 210, 17, 6],
    "Rn": ["Rn", "Radon", 222, 18, 6],
    # Period 7 - Including Actinides
    "Fr": ["Fr", "Francium", 223, 1, 7], "Ra": ["Ra", "Radium", 226, 2, 7],
    "Ac": ["Ac", "Actinium", 227, 3, 7], "Th": ["Th", "Thorium", 232.04, 101, 7],
    "Pa": ["Pa", "Protactinium", 231.04, 101, 7], "U": ["U", "Uranium", 238.03, 101, 7],
    "Np": ["Np", "Neptunium", 237, 101, 7], "Pu": ["Pu", "Plutonium", 244, 101, 7],
    "Am": ["Am", "Americium", 243, 101, 7], "Cm": ["Cm", "Curium", 247, 101, 7],
    "Bk": ["Bk", "Berkelium", 247, 101, 7], "Cf": ["Cf", "Californium", 251, 101, 7],
    "Es": ["Es", "Einsteinium", 252, 101, 7], "Fm": ["Fm", "Fermium", 257, 101, 7],
    "Md": ["Md", "Mendelevium", 258, 101, 7], "No": ["No", "Nobelium", 259, 101, 7],
    "Lr": ["Lr", "Lawrencium", 262, 3, 7],
    "Rf": ["Rf", "Rutherfordium", 267, 4, 7], "Db": ["Db", "Dubnium", 268, 5, 7],
    "Sg": ["Sg", "Seaborgium", 269, 6, 7], "Bh": ["Bh", "Bohrium", 270, 7, 7],
    "Hs": ["Hs", "Hassium", 269, 8, 7], "Mt": ["Mt", "Meitnerium", 278, 9, 7],
    "Ds": ["Ds", "Darmstadtium", 281, 10, 7], "Rg": ["Rg", "Roentgenium", 282, 11, 7],
    "Cn": ["Cn", "Copernicium", 285, 12, 7], "Nh": ["Nh", "Nihonium", 286, 13, 7],
    "Fl": ["Fl", "Flerovium", 289, 14, 7], "Mc": ["Mc", "Moscovium", 290, 15, 7],
    "Lv": ["Lv", "Livermorium", 293, 16, 7], "Ts": ["Ts", "Tennessine", 294, 17, 7],
    "Og": ["Og", "Oganesson", 294, 18, 7]
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
    "Ethanol": {"formula": "C2H5OH", "mass": 46.07},
    "Aspirin": {"formula": "C9H8O4", "mass": 180.157},
    "Caffeine": {"formula": "C8H10N4O2", "mass": 194.19}
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
                "description": "Perfect tetrahedral molecule with Td symmetry."},
    "ammonia": {"name": "Ammonia", "formula": "NH₃", "point_group": "C₃v",
                "elements": [{"type": "C₃", "axis": [0, 1, 0], "description": "3-fold rotation axis"}],
                "description": "Trigonal pyramidal with C₃v symmetry."}
}

# Molecule 3D structures
MOLECULE_STRUCTURES = {
    "water": "O", "benzene": "c1ccccc1", "methane": "C",
    "ammonia": "N", "carbon_dioxide": "O=C=O", "ethanol": "CCO"
}

# ========== UTILITY FUNCTIONS ==========
def get_theme_colors():
    """Get theme-adaptive colors based on Streamlit theme"""
    try:
        is_dark = st.get_option("theme.base") == "dark"
    except:
        is_dark = False
    
    if is_dark:
        return {
            "bg_color": "#1e1e2f",
            "card_bg": "#2d2d3f",
            "text_color": "#ffffff",
            "text_secondary": "#b0b0b0",
            "border_color": "#404060",
            "result_bg": "#2d2d3f",
            "result_border": "#667eea",
            "success_bg": "#1a3a2a",
            "success_text": "#00ff88",
            "warning_bg": "#3a2a1a",
            "warning_text": "#ffaa44",
            "info_bg": "#1a2a3a",
            "info_text": "#44aaff"
        }
    else:
        return {
            "bg_color": "#ffffff",
            "card_bg": "#f8f9fa",
            "text_color": "#1e1e2f",
            "text_secondary": "#666666",
            "border_color": "#e0e0e0",
            "result_bg": "#e8f4f8",
            "result_border": "#667eea",
            "success_bg": "#d4edda",
            "success_text": "#155724",
            "warning_bg": "#fff3cd",
            "warning_text": "#856404",
            "info_bg": "#d1ecf1",
            "info_text": "#0c5460"
        }

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
        "n2+h2->nh3": "N₂ + 3H₂ → 2NH₃",
        "fe+o2->fe2o3": "4Fe + 3O₂ → 2Fe₂O₃",
        "c2h5oh+o2->co2+h2o": "C₂H₅OH + 3O₂ → 2CO₂ + 3H₂O"
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

# ========== DISPLAY RESULT CARD WITH ADAPTIVE THEMING ==========
def display_result_card(mass, molar_mass, moles, molecules):
    """Display result card with adaptive theming"""
    colors = get_theme_colors()
    
    # Inject CSS variables for theming
    st.markdown(f"""
    <style>
        .result-card {{
            --result-bg: {colors['result_bg']};
            --result-border: {colors['result_border']};
            --text-color: {colors['text_color']};
        }}
    </style>
    """, unsafe_allow_html=True)
    
    result_html = f"""
    <div class="result-card">
        <strong>📊 Result:</strong><br>
        Mass: {mass:.4f} g<br>
        Molar Mass: {molar_mass:.4f} g/mol<br>
        <strong>Moles = {moles:.6f} mol</strong><br>
        Molecules: {molecules:.2e}
    </div>
    """
    st.markdown(result_html, unsafe_allow_html=True)

# ========== MAIN APPLICATION ==========
def main():
    # Get theme colors
    colors = get_theme_colors()
    
    # Inject theme variables
    st.markdown(f"""
    <style>
        :root {{
            --result-bg: {colors['result_bg']};
            --result-border: {colors['result_border']};
            --card-bg: {colors['card_bg']};
            --text-color: {colors['text_color']};
            --border-color: {colors['border_color']};
        }}
    </style>
    """, unsafe_allow_html=True)
    
    # Header
    st.markdown("""
    <div class="main-header">
        <h1 style="color: white; margin: 0;">🧪 Ojo Chem Pro</h1>
        <p style="color: white; margin: 0; opacity: 0.9;">Advanced 3D Chemistry Lab with Mass-Mole Converter</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Emergency banner
    with st.expander("🚨 Emergency Information", expanded=False):
        st.warning("**Poison Control:** 1-800-222-1222 | **CHEMTREC:** 1-800-424-9300\n\n**First Aid:** Flush with water for 15 minutes, seek immediate medical attention.")
    
    # Sidebar
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
        compound_option = st.radio("Compound:", ["Common", "Custom Formula", "Periodic Table"], horizontal=True)
        
        molar_mass = None
        formula_display = ""
        
        if compound_option == "Common":
            compound_name = st.selectbox("Select:", list(COMMON_COMPOUNDS.keys()))
            formula_display = COMMON_COMPOUNDS[compound_name]["formula"]
            molar_mass = COMMON_COMPOUNDS[compound_name]["mass"]
            st.info(f"**Formula:** {format_formula(formula_display)} | **Molar Mass:** {molar_mass:.3f} g/mol")
        
        elif compound_option == "Custom Formula":
            formula_input = st.text_input("Formula:", placeholder="e.g., H2O, C6H12O6, KMnO4")
            if formula_input:
                mass_calc, _ = calculate_molar_mass(formula_input.upper())
                if mass_calc:
                    molar_mass = mass_calc
                    formula_display = formula_input
                    st.success(f"**Molar Mass:** {molar_mass:.3f} g/mol")
        
        else:  # Periodic Table - Build compound from elements
            st.markdown("#### Build from Elements")
            
            # Element selection
            element_list = sorted([(symbol, data[1], data[2]) for symbol, data in PERIODIC_TABLE.items()])
            element_options = [f"{symbol} - {name} ({atomic_mass:.3f} g/mol)" for symbol, name, atomic_mass in element_list]
            
            selected_element = st.selectbox("Select element:", element_options)
            selected_symbol = selected_element.split(" - ")[0]
            
            col1, col2 = st.columns(2)
            with col1:
                element_count = st.number_input("Count:", min_value=1, max_value=100, value=1, step=1)
            with col2:
                if st.button("➕ Add", use_container_width=True):
                    if "compound_elements" not in st.session_state:
                        st.session_state.compound_elements = {}
                    st.session_state.compound_elements[selected_symbol] = st.session_state.compound_elements.get(selected_symbol, 0) + element_count
                    st.rerun()
            
            # Display current composition
            if "compound_elements" in st.session_state and st.session_state.compound_elements:
                st.markdown("**Current composition:**")
                for element, count in st.session_state.compound_elements.items():
                    st.markdown(f"- {element}: {count}")
                
                col1, col2 = st.columns(2)
                with col1:
                    if st.button("🗑️ Clear", use_container_width=True):
                        st.session_state.compound_elements = {}
                        st.rerun()
                with col2:
                    if st.button("🔨 Calculate", use_container_width=True):
                        total_mass = 0
                        formula_parts = []
                        for element, count in st.session_state.compound_elements.items():
                            if element in ATOMIC_MASSES:
                                total_mass += ATOMIC_MASSES[element] * count
                                formula_parts.append(f"{element}{count if count > 1 else ''}")
                        if total_mass > 0:
                            molar_mass = total_mass
                            formula_display = "".join(formula_parts)
                            st.success(f"**Formula:** {format_formula(formula_display)}\n**Molar Mass:** {molar_mass:.3f} g/mol")
        
        # Conversion calculator with adaptive result card
        if molar_mass and molar_mass > 0:
            st.divider()
            if conversion_type == "Mass → Moles":
                mass_input = st.number_input("Mass (grams):", min_value=0.0, value=1.0, step=0.1, format="%.4f")
                if mass_input > 0:
                    moles = mass_to_moles(mass_input, molar_mass)
                    molecules = moles * 6.022e23
                    display_result_card(mass_input, molar_mass, moles, molecules)
            else:
                moles_input = st.number_input("Moles:", min_value=0.0, value=1.0, step=0.1, format="%.4f")
                if moles_input > 0:
                    mass = moles_to_mass(moles_input, molar_mass)
                    molecules = moles_input * 6.022e23
                    st.markdown(f"""
                    <div class="result-card">
                        <strong>📊 Result:</strong><br>
                        Moles: {moles_input:.4f} mol<br>
                        Molar Mass: {molar_mass:.4f} g/mol<br>
                        <strong>Mass = {mass:.4f} g</strong><br>
                        Molecules: {molecules:.2e}
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
    
    # Main content area
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
    
    # RIGHT COLUMN - Chat interface
    with col_right:
        st.markdown("### 💬 Chemistry Assistant")
        
        if "chat_messages" not in st.session_state:
            st.session_state.chat_messages = [
                {"role": "assistant", "content": "👋 Hi! Ask me about molecules, symmetry, or chemistry!"}
            ]
        
        # Display chat messages
        chat_container = st.container()
        with chat_container:
            st.markdown(f'<div style="height: 450px; overflow-y: auto; padding: 1rem; background: {colors["card_bg"]}; border-radius: 10px; margin-bottom: 1rem; border: 1px solid {colors["border_color"]};">', unsafe_allow_html=True)
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
    
    # Chat input at bottom
    st.divider()
    user_input = st.chat_input("Ask me about chemistry...", key="chemistry_chat_input_bottom")
    
    if user_input:
        st.session_state.chat_messages.append({"role": "user", "content": user_input})
        with st.spinner("Thinking..."):
            response = get_chatbot_response(user_input, safe_mode)
        st.session_state.chat_messages.append({"role": "assistant", "content": response})
        st.rerun()

if __name__ == "__main__":
    main()
