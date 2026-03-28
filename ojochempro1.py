"""
Ojo Chem Pro - Advanced 3D Chemistry Lab with Symmetry Analysis
Streamlit Deployment Version - Compatible with Python 3.9+
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import requests
import re
from datetime import datetime
import json

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
    
    .chat-message {
        padding: 0.75rem;
        border-radius: 10px;
        margin-bottom: 0.5rem;
        max-width: 80%;
    }
    
    .user-message {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        float: right;
        clear: both;
    }
    
    .bot-message {
        background: #f0f2f6;
        color: #1e1e2f;
        float: left;
        clear: both;
        border: 1px solid #e0e0e0;
    }
    
    .symmetry-badge {
        display: inline-block;
        padding: 0.2rem 0.6rem;
        margin: 0.2rem;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border-radius: 20px;
        font-size: 0.8rem;
        cursor: pointer;
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
    
    /* 3D container styling */
    .mol-container {
        width: 100%;
        height: 500px;
        position: relative;
        background: #1a1a2e;
        border-radius: 10px;
        overflow: hidden;
    }
</style>
""", unsafe_allow_html=True)

# ========== CHEMICAL DATABASE ==========
ATOMIC_MASSES = {
    'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81, 'C': 12.011,
    'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.18, 'Na': 22.99, 'Mg': 24.305,
    'Al': 26.982, 'Si': 28.086, 'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.95,
    'K': 39.098, 'Ca': 40.078, 'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996,
    'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.63, 'As': 74.922, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798,
    'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224, 'Nb': 92.906, 'Mo': 95.95,
    'Tc': 98, 'Ru': 101.07, 'Rh': 102.91, 'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41,
    'In': 114.82, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.6, 'I': 126.90, 'Xe': 131.29,
    'Cs': 132.91, 'Ba': 137.33, 'La': 138.91, 'Ce': 140.12, 'Pr': 140.91, 'Nd': 144.24,
    'Pm': 145, 'Sm': 150.36, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93, 'Dy': 162.5,
    'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93, 'Yb': 173.05, 'Lu': 174.97, 'Hf': 178.49,
    'Ta': 180.95, 'W': 183.84, 'Re': 186.21, 'Os': 190.23, 'Ir': 192.22, 'Pt': 195.08,
    'Au': 196.97, 'Hg': 200.59, 'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.98, 'Th': 232.04,
    'Pa': 231.04, 'U': 238.03
}

# Molecule symmetry database
MOLECULE_SYMMETRY = {
    "water": {
        "name": "Water",
        "formula": "H₂O",
        "point_group": "C₂v",
        "elements": [
            {"type": "C₂", "axis": [0, 1, 0], "description": "2-fold rotation axis through O atom"},
            {"type": "σᵥ", "plane": [1, 0, 0], "description": "Vertical mirror plane (molecular plane)"},
            {"type": "σᵥ'", "plane": [0, 0, 1], "description": "Vertical mirror plane (perpendicular)"}
        ],
        "description": "Bent molecule with C₂v symmetry. Has a C₂ axis through the oxygen atom and two mirror planes."
    },
    "benzene": {
        "name": "Benzene",
        "formula": "C₆H₆",
        "point_group": "D₆h",
        "elements": [
            {"type": "C₆", "axis": [0, 1, 0], "description": "6-fold principal axis through ring center"},
            {"type": "C₂", "axis": [1, 0, 0], "description": "2-fold axes through opposite carbons"},
            {"type": "σₕ", "plane": [0, 1, 0], "description": "Horizontal mirror plane (molecular plane)"},
            {"type": "σᵥ", "plane": [1, 0, 0], "description": "Vertical mirror planes through atoms"},
            {"type": "i", "center": [0, 0, 0], "description": "Inversion center at ring center"}
        ],
        "description": "Aromatic ring with D₆h symmetry. High symmetry with multiple rotation axes and mirror planes."
    },
    "methane": {
        "name": "Methane",
        "formula": "CH₄",
        "point_group": "Td",
        "elements": [
            {"type": "C₃", "axis": [1, 1, 1], "description": "3-fold axes through C-H bonds"},
            {"type": "C₂", "axis": [1, 0, 0], "description": "2-fold axes through midpoints of edges"},
            {"type": "σd", "plane": [1, 1, 0], "description": "Dihedral mirror planes"}
        ],
        "description": "Perfect tetrahedral molecule with Td symmetry. Highest symmetry among simple molecules."
    },
    "ammonia": {
        "name": "Ammonia",
        "formula": "NH₃",
        "point_group": "C₃v",
        "elements": [
            {"type": "C₃", "axis": [0, 1, 0], "description": "3-fold rotation axis through N atom"},
            {"type": "σᵥ", "plane": [1, 0, 0], "description": "Vertical mirror planes through N and H"}
        ],
        "description": "Trigonal pyramidal with C₃v symmetry. Has a C₃ axis through the nitrogen atom."
    },
    "carbon_dioxide": {
        "name": "Carbon Dioxide",
        "formula": "CO₂",
        "point_group": "D∞h",
        "elements": [
            {"type": "C∞", "axis": [0, 1, 0], "description": "Infinite rotation axis along molecule"},
            {"type": "C₂", "axis": [1, 0, 0], "description": "2-fold axes perpendicular to principal axis"},
            {"type": "σₕ", "plane": [0, 1, 0], "description": "Horizontal mirror plane (perpendicular to axis)"},
            {"type": "σᵥ", "plane": [1, 0, 0], "description": "Vertical mirror planes containing axis"},
            {"type": "i", "center": [0, 0, 0], "description": "Inversion center at carbon atom"}
        ],
        "description": "Linear molecule with D∞h symmetry. Has infinite rotation axis and inversion center."
    },
    "ethanol": {
        "name": "Ethanol",
        "formula": "C₂H₅OH",
        "point_group": "C₁",
        "elements": [],
        "description": "Asymmetric molecule with only identity symmetry (C₁). No symmetry elements besides identity."
    }
}

# Molecule 3D structure data (3D coordinates for visualization)
MOLECULE_COORDINATES = {
    "water": {
        "atoms": [
            {"element": "O", "x": 0, "y": 0, "z": 0, "color": "#ff6666", "size": 0.5},
            {"element": "H", "x": 0.8, "y": 0.6, "z": 0, "color": "#cccccc", "size": 0.3},
            {"element": "H", "x": 0.8, "y": -0.6, "z": 0, "color": "#cccccc", "size": 0.3}
        ],
        "bonds": [[0, 1], [0, 2]]
    },
    "methane": {
        "atoms": [
            {"element": "C", "x": 0, "y": 0, "z": 0, "color": "#888888", "size": 0.5},
            {"element": "H", "x": 0.8, "y": 0.8, "z": 0.8, "color": "#cccccc", "size": 0.3},
            {"element": "H", "x": 0.8, "y": -0.8, "z": -0.8, "color": "#cccccc", "size": 0.3},
            {"element": "H", "x": -0.8, "y": 0.8, "z": -0.8, "color": "#cccccc", "size": 0.3},
            {"element": "H", "x": -0.8, "y": -0.8, "z": 0.8, "color": "#cccccc", "size": 0.3}
        ],
        "bonds": [[0, 1], [0, 2], [0, 3], [0, 4]]
    },
    "ammonia": {
        "atoms": [
            {"element": "N", "x": 0, "y": 0, "z": 0, "color": "#44aaff", "size": 0.5},
            {"element": "H", "x": 0.8, "y": 0.8, "z": 0, "color": "#cccccc", "size": 0.3},
            {"element": "H", "x": 0.8, "y": -0.4, "z": 0.7, "color": "#cccccc", "size": 0.3},
            {"element": "H", "x": 0.8, "y": -0.4, "z": -0.7, "color": "#cccccc", "size": 0.3}
        ],
        "bonds": [[0, 1], [0, 2], [0, 3]]
    },
    "benzene": {
        "atoms": [],
        "bonds": []
    }
}

# Generate benzene coordinates
def generate_benzene_coordinates():
    atoms = []
    bonds = []
    radius = 1.2
    
    for i in range(6):
        angle = (i * 60) * np.pi / 180
        x = np.cos(angle) * radius
        z = np.sin(angle) * radius
        atoms.append({"element": "C", "x": x, "y": 0, "z": z, "color": "#888888", "size": 0.4})
        atoms.append({"element": "H", "x": x * 1.6, "y": 0, "z": z * 1.6, "color": "#cccccc", "size": 0.28})
        bonds.append([i*2, i*2 + 1])
    
    for i in range(6):
        bonds.append([i*2, ((i+1)%6)*2])
    
    return {"atoms": atoms, "bonds": bonds}

MOLECULE_COORDINATES["benzene"] = generate_benzene_coordinates()

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
            # Get element symbol
            if i + 1 < len(formula) and formula[i + 1].islower():
                element = formula[i:i + 2]
                i += 2
            else:
                element = formula[i]
                i += 1
            
            # Get count
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

def balance_equation(equation):
    """Balance chemical equations"""
    balanced_equations = {
        "h2+o2->h2o": "2H₂ + O₂ → 2H₂O",
        "ch4+o2->co2+h2o": "CH₄ + 2O₂ → CO₂ + 2H₂O",
        "c6h12o6+o2->co2+h2o": "C₆H₁₂O₆ + 6O₂ → 6CO₂ + 6H₂O",
        "na+cl2->nacl": "2Na + Cl₂ → 2NaCl",
        "n2+h2->nh3": "N₂ + 3H₂ → 2NH₃",
        "fe+o2->fe2o3": "4Fe + 3O₂ → 2Fe₂O₃",
        "c2h5oh+o2->co2+h2o": "C₂H₅OH + 3O₂ → 2CO₂ + 3H₂O",
        "kmno4+hcl->kcl+mncl2+cl2+h2o": "2KMnO₄ + 16HCl → 2KCl + 2MnCl₂ + 5Cl₂ + 8H₂O"
    }
    
    key = equation.lower().replace(" ", "").replace("->", "->")
    if key in balanced_equations:
        return balanced_equations[key]
    return None

def search_pubchem(compound):
    """Search PubChem API for compound information"""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound}/property/MolecularFormula,MolecularWeight,InChIKey/JSON"
        response = requests.get(url, timeout=5)
        
        if response.status_code == 200:
            data = response.json()
            if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
                prop = data["PropertyTable"]["Properties"][0]
                return {
                    "formula": prop.get("MolecularFormula", "N/A"),
                    "weight": prop.get("MolecularWeight", "N/A"),
                    "inchikey": prop.get("InChIKey", "N/A")
                }
    except:
        pass
    return None

# ========== 3D MOLECULE VISUALIZATION WITH PLOTLY ==========
def create_3d_molecule_plotly(molecule_name):
    """Create 3D molecule visualization using Plotly"""
    molecule_name = molecule_name.lower()
    
    if molecule_name not in MOLECULE_COORDINATES:
        return None
    
    coords = MOLECULE_COORDINATES[molecule_name]
    
    fig = go.Figure()
    
    # Add atoms as spheres
    for atom in coords["atoms"]:
        fig.add_trace(go.Scatter3d(
            x=[atom["x"]],
            y=[atom["y"]],
            z=[atom["z"]],
            mode='markers+text',
            marker=dict(
                size=atom["size"] * 20,
                color=atom["color"],
                symbol='circle',
                line=dict(color='white', width=1)
            ),
            text=[atom["element"]],
            textposition='top center',
            name=f"{atom['element']}",
            hoverinfo='text',
            hovertext=f"{atom['element']}<br>Position: ({atom['x']:.2f}, {atom['y']:.2f}, {atom['z']:.2f})"
        ))
    
    # Add bonds as lines
    for bond in coords["bonds"]:
        atom1 = coords["atoms"][bond[0]]
        atom2 = coords["atoms"][bond[1]]
        fig.add_trace(go.Scatter3d(
            x=[atom1["x"], atom2["x"]],
            y=[atom1["y"], atom2["y"]],
            z=[atom1["z"], atom2["z"]],
            mode='lines',
            line=dict(color='#ccccaa', width=4),
            showlegend=False,
            hoverinfo='none'
        ))
    
    # Update layout
    fig.update_layout(
        title=f"{MOLECULE_SYMMETRY.get(molecule_name, {}).get('name', molecule_name.title())} - 3D Structure",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z",
            bgcolor='#1a1a2e',
            xaxis=dict(color='white', gridcolor='#335588'),
            yaxis=dict(color='white', gridcolor='#335588'),
            zaxis=dict(color='white', gridcolor='#335588'),
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5)
            )
        ),
        paper_bgcolor='#1a1a2e',
        plot_bgcolor='#1a1a2e',
        font=dict(color='white'),
        height=500,
        hovermode='closest'
    )
    
    return fig

# ========== SYMMETRY VISUALIZATION ==========
def create_symmetry_visualization(molecule_name):
    """Create symmetry element visualization using Plotly"""
    sym_data = MOLECULE_SYMMETRY.get(molecule_name.lower())
    if not sym_data or not sym_data.get("elements"):
        return None
    
    fig = go.Figure()
    
    # Add sphere for molecule center
    fig.add_trace(go.Scatter3d(
        x=[0], y=[0], z=[0],
        mode='markers',
        marker=dict(size=6, color='white', symbol='circle'),
        name='Molecule Center',
        hoverinfo='text',
        hovertext='Center of mass'
    ))
    
    # Add symmetry elements
    colors = {'C₂': '#ff6666', 'C₃': '#ffaa66', 'C₆': '#ffff66', 
              'σᵥ': '#66ff66', 'σₕ': '#66ff66', 'σd': '#66ff66', 
              'i': '#ffaa44', 'C∞': '#ff8888'}
    
    for element in sym_data["elements"]:
        if "axis" in element["type"] or "C" in element["type"]:
            # Draw rotation axis
            axis = element.get("axis", [0, 1, 0])
            length = 2.5
            fig.add_trace(go.Scatter3d(
                x=[-axis[0]*length, axis[0]*length],
                y=[-axis[1]*length, axis[1]*length],
                z=[-axis[2]*length, axis[2]*length],
                mode='lines',
                line=dict(color=colors.get(element["type"], '#ff6666'), width=6),
                name=f"{element['type']} Axis",
                hoverinfo='text',
                hovertext=element["description"]
            ))
            
            # Add rotation indicator (circle)
            if "C" in element["type"] and element["type"] != "C∞":
                order = int(re.search(r'\d+', element["type"]).group()) if re.search(r'\d+', element["type"]) else 2
                radius = 1.0
                theta = np.linspace(0, 2*np.pi, 100)
                
                # Find perpendicular plane
                perp_axis = [1, 0, 0] if abs(axis[1]) > 0.8 else [0, 0, 1]
                perp_axis = np.array(perp_axis) - np.dot(perp_axis, axis) * np.array(axis)
                perp_axis = perp_axis / np.linalg.norm(perp_axis)
                
                x_circle = radius * np.cos(theta) * perp_axis[0]
                y_circle = radius * np.cos(theta) * perp_axis[1]
                z_circle = radius * np.sin(theta) * perp_axis[2] if perp_axis[2] != 0 else radius * np.cos(theta) * perp_axis[1]
                
                fig.add_trace(go.Scatter3d(
                    x=x_circle, y=y_circle, z=z_circle,
                    mode='lines',
                    line=dict(color='#ff8888', width=2, dash='dash'),
                    showlegend=False,
                    hoverinfo='none'
                ))
                
        elif "plane" in element["type"] or "σ" in element["type"]:
            # Draw mirror plane (semi-transparent surface)
            plane = element.get("plane", [1, 0, 0])
            size = 2.0
            
            # Create grid for plane
            u, v = np.mgrid[-size:size:20j, -size:size:20j]
            if plane[0] != 0:
                # Plane perpendicular to x-axis
                x = np.zeros_like(u)
                y = u
                z = v
            elif plane[2] != 0:
                # Plane perpendicular to z-axis
                x = u
                y = v
                z = np.zeros_like(u)
            else:
                # Plane perpendicular to y-axis
                x = u
                y = np.zeros_like(u)
                z = v
            
            fig.add_trace(go.Surface(
                x=x, y=y, z=z,
                opacity=0.25,
                colorscale=[[0, 'green'], [1, 'lightgreen']],
                showscale=False,
                name=f"{element['type']} Plane",
                hoverinfo='text',
                hovertext=element["description"]
            ))
            
        elif "center" in element["type"] or "i" in element["type"]:
            # Draw inversion center
            fig.add_trace(go.Scatter3d(
                x=[0], y=[0], z=[0],
                mode='markers',
                marker=dict(size=12, color='#ffaa44', symbol='circle-open', 
                           line=dict(color='#ffaa44', width=3)),
                name='Inversion Center (i)',
                hoverinfo='text',
                hovertext=element["description"]
            ))
    
    # Update layout
    fig.update_layout(
        title=f"Symmetry Elements of {sym_data['name']} (Point Group: {sym_data['point_group']})",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z",
            bgcolor='#1a1a2e',
            xaxis=dict(color='white', gridcolor='#335588'),
            yaxis=dict(color='white', gridcolor='#335588'),
            zaxis=dict(color='white', gridcolor='#335588'),
            camera=dict(
                eye=dict(x=1.8, y=1.5, z=1.2)
            )
        ),
        paper_bgcolor='#1a1a2e',
        plot_bgcolor='#1a1a2e',
        font=dict(color='white'),
        height=500
    )
    
    return fig

# ========== CHATBOT RESPONSE ENGINE ==========
def get_chatbot_response(query, safe_mode=True):
    """Generate response for user queries"""
    query_lower = query.lower()
    
    # Safety check
    hazardous_keywords = ['explosive', 'cyanide', 'mustard gas', 'weapon', 'toxic gas', 'nerve agent']
    if safe_mode and any(keyword in query_lower for keyword in hazardous_keywords):
        return "⚠️ **Safety Notice**: This topic involves hazardous materials and is blocked in Safe Mode. Please consult with a teacher or guardian."
    
    # Molecule visualization
    for molecule in MOLECULE_SYMMETRY.keys():
        if molecule in query_lower:
            sym = MOLECULE_SYMMETRY[molecule]
            return f"🔬 **{sym['name']}**\n\n" \
                   f"**Formula:** {sym['formula']}\n" \
                   f"**Point Group:** {sym['point_group']}\n\n" \
                   f"**Description:** {sym['description']}\n\n" \
                   f"**Symmetry Elements:**\n" + \
                   "\n".join([f"• {elem['description']}" for elem in sym['elements']]) + \
                   f"\n\n💡 **Tip:** Check the 3D Viewer tab to see the molecule structure, and the Symmetry Analysis tab for 3D visualization of symmetry elements!"
    
    # Equation balancing
    if "balance" in query_lower or "->" in query_lower:
        balanced = balance_equation(query)
        if balanced:
            return f"✅ **Balanced Equation:**\n\n{balanced}\n\n📝 **Steps:**\n1. Count atoms on each side\n2. Add coefficients to balance\n3. Verify all atoms are equal"
        return "🤔 I can help balance equations! Try examples like:\n• H2 + O2 -> H2O\n• CH4 + O2 -> CO2 + H2O\n• C6H12O6 + O2 -> CO2 + H2O\n• KMnO4 + HCl -> KCl + MnCl2 + Cl2 + H2O"
    
    # Molar mass calculation
    if "molar mass" in query_lower or "molecular weight" in query_lower:
        # Extract formula
        formula_match = re.search(r'([A-Z][a-z]?\d*)+', query.upper())
        if formula_match:
            formula = formula_match.group()
            mass, elements = calculate_molar_mass(formula)
            if mass:
                return f"📊 **Molar Mass of {format_formula(formula)}**\n\n**Molar Mass:** {mass:.4f} g/mol\n\n**Composition:**\n" + "\n".join(elements[:5])
        return "Please provide a chemical formula (e.g., 'Molar mass of H2O' or 'C6H12O6')"
    
    # PubChem search
    if "pubchem" in query_lower or "search for" in query_lower:
        compound = query_lower.replace("pubchem", "").replace("search for", "").strip()
        if compound:
            result = search_pubchem(compound)
            if result:
                return f"🔍 **PubChem Results for {compound.title()}**\n\n" \
                       f"📝 **Formula:** {result['formula']}\n" \
                       f"⚖️ **Molecular Weight:** {result['weight']} g/mol\n" \
                       f"🔑 **InChIKey:** {result['inchikey']}\n\n" \
                       f"[View on PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/{compound})"
            else:
                return f"🔍 Compound '{compound}' not found. Try checking the spelling or use a different name."
        return "Please specify a compound to search (e.g., 'Search PubChem for aspirin')"
    
    # Point group query
    if "point group" in query_lower or "symmetry" in query_lower:
        for molecule, data in MOLECULE_SYMMETRY.items():
            if molecule in query_lower:
                return f"**{data['name']}** belongs to point group **{data['point_group']}**.\n\n{data['description']}\n\n**Symmetry elements:**\n" + "\n".join([f"• {elem['description']}" for elem in data['elements']])
    
    # General chemistry help
    return f"🧪 **Chemistry Assistant**\n\nI can help you with:\n\n" \
           f"• **3D Molecule Viewer** - View molecules in 3D (use the 3D Viewer tab)\n" \
           f"• **Symmetry Analysis** - Visualize symmetry elements (use Symmetry Analysis tab)\n" \
           f"• **Equation Balancing** - Balance chemical equations\n" \
           f"• **Molar Mass Calculator** - Calculate molecular weights\n" \
           f"• **PubChem Integration** - Search chemical compounds\n\n" \
           f"**Try these commands:**\n" \
           f"• 'Show me benzene'\n" \
           f"• 'What is the symmetry of water?'\n" \
           f"• 'Balance: H2 + O2 -> H2O'\n" \
           f"• 'Molar mass of glucose'\n" \
           f"• 'Search PubChem for aspirin'\n\n" \
           f"**Current Mode:** {'🔒 Safe Mode (Educational)' if safe_mode else '🔬 Expert Mode'}"
    
# ========== MAIN APPLICATION ==========
def main():
    # Header
    st.markdown("""
    <div class="main-header">
        <h1 style="color: white; margin: 0;">🧪 Ojo Chem Pro</h1>
        <p style="color: white; margin: 0; opacity: 0.9;">Advanced 3D Chemistry Lab with Symmetry Analysis</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Sidebar for mode and controls
    with st.sidebar:
        st.markdown("### 🎮 Controls")
        
        # Mode toggle
        safe_mode = st.toggle("🔒 Safe Mode", value=True, 
                              help="Safe Mode: Educational content only. Expert Mode: Advanced features with safety warnings.")
        
        if not safe_mode:
            st.warning("⚠️ Expert Mode: Advanced content with safety warnings. Adult supervision recommended.")
        
        st.divider()
        
        # Quick molecule library
        st.markdown("### 📚 Quick Molecules")
        
        molecule_buttons = {
            "💧 Water": "water",
            "🔄 Benzene": "benzene",
            "🔥 Methane": "methane",
            "🧪 Ammonia": "ammonia",
            "💨 CO₂": "carbon_dioxide",
            "🍷 Ethanol": "ethanol"
        }
        
        cols = st.columns(3)
        for i, (label, mol_key) in enumerate(molecule_buttons.items()):
            with cols[i % 3]:
                if st.button(label, key=f"quick_{mol_key}", use_container_width=True):
                    st.session_state['selected_molecule'] = mol_key
                    st.session_state['active_tab'] = "Viewer"
                    st.rerun()
        
        st.divider()
        
        # Emergency information (minimized)
        with st.expander("🚨 Emergency Info", expanded=False):
            st.warning("""
            **Poison Control:** 1-800-222-1222  
            **CHEMTREC:** 1-800-424-9300  
            
            **First Aid:** Flush with water for 15 minutes, seek immediate medical attention.
            """)
        
        st.divider()
        
        # About
        st.markdown("### ℹ️ About")
        st.info("""
        **Ojo Chem Pro** helps you explore:
        - 3D molecular structures
        - Symmetry elements and point groups
        - Chemical calculations
        - Equation balancing
        """)
        
        # Stats
        st.markdown(f"**Molecules:** {len(MOLECULE_SYMMETRY)}")
        st.markdown(f"**Elements:** {len(ATOMIC_MASSES)}")
    
    # Main content area with tabs
    tab1, tab2, tab3, tab4 = st.tabs(["🔬 3D Viewer", "🔄 Symmetry Analysis", "⚙️ Tools", "💬 Chemistry Chat"])
    
    # Tab 1: 3D Viewer
    with tab1:
        st.markdown("#### 3D Molecular Structure")
        
        # Molecule selector
        molecule_options = list(MOLECULE_COORDINATES.keys())
        default_index = 0
        if 'selected_molecule' in st.session_state and st.session_state['selected_molecule'] in molecule_options:
            default_index = molecule_options.index(st.session_state['selected_molecule'])
        
        molecule_name = st.selectbox(
            "Select Molecule",
            options=molecule_options,
            format_func=lambda x: MOLECULE_SYMMETRY.get(x, {}).get('name', x.title()),
            index=default_index,
            key="viewer_selector"
        )
        
        # Create and display 3D visualization
        fig = create_3d_molecule_plotly(molecule_name)
        if fig:
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("3D visualization loading...")
        
        # Molecule information
        if molecule_name in MOLECULE_SYMMETRY:
            sym = MOLECULE_SYMMETRY[molecule_name]
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Formula", sym['formula'])
                st.metric("Point Group", sym['point_group'])
            with col2:
                st.metric("Atoms", len(MOLECULE_COORDINATES.get(molecule_name, {}).get('atoms', [])))
                st.metric("Bonds", len(MOLECULE_COORDINATES.get(molecule_name, {}).get('bonds', [])))
            
            st.info(sym['description'])
    
    # Tab 2: Symmetry Analysis
    with tab2:
        st.markdown("#### 🔄 Symmetry Elements Visualization")
        
        sym_molecule = st.selectbox(
            "Select Molecule for Symmetry Analysis",
            options=list(MOLECULE_SYMMETRY.keys()),
            format_func=lambda x: MOLECULE_SYMMETRY[x]['name'],
            key="sym_selector"
        )
        
        # Display symmetry visualization
        fig = create_symmetry_visualization(sym_molecule)
        if fig:
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("This molecule has no symmetry elements (C₁ point group) or visualization is being generated.")
        
        # Symmetry information
        sym_data = MOLECULE_SYMMETRY.get(sym_molecule, {})
        st.markdown(f"### {sym_data.get('name', 'N/A')}")
        st.markdown(f"**Formula:** {sym_data.get('formula', 'N/A')}")
        st.markdown(f"**Point Group:** `{sym_data.get('point_group', 'N/A')}`")
        
        if sym_data.get('elements'):
            st.markdown("#### Symmetry Elements")
            for elem in sym_data['elements']:
                st.markdown(f"- **{elem['type']}**: {elem['description']}")
        
        with st.expander("📖 Understanding Symmetry"):
            st.markdown("""
            **Rotation Axes (Cₙ):** Axis around which molecule rotates to appear identical.
            
            **Mirror Planes (σ):** Planes that reflect molecule onto itself.
            - **σᵥ** - Vertical mirror plane (contains principal axis)
            - **σₕ** - Horizontal mirror plane (perpendicular to principal axis)
            - **σd** - Dihedral mirror plane (bisects C₂ axes)
            
            **Inversion Center (i):** Point where (x,y,z) → (-x,-y,-z).
            
            **Point Groups:** Classification of molecular symmetry. Common point groups:
            - **C₁** - No symmetry
            - **C₂v** - Water, SO₂
            - **C₃v** - Ammonia, CHCl₃
            - **Td** - Methane, CCl₄
            - **D₆h** - Benzene
            - **D∞h** - CO₂, H₂
            """)
    
    # Tab 3: Tools
    with tab3:
        st.markdown("#### ⚙️ Chemistry Tools")
        
        tool_tabs = st.tabs(["📊 Molar Mass", "⚖️ Equation Balancer", "🔍 PubChem Search"])
        
        # Molar Mass Calculator
        with tool_tabs[0]:
            st.markdown("##### Calculate Molar Mass")
            formula_input = st.text_input("Chemical Formula", placeholder="e.g., H2O, C6H12O6, NaCl", key="molar_formula")
            
            if st.button("Calculate", key="calc_molar"):
                if formula_input:
                    with st.spinner("Calculating..."):
                        mass, elements = calculate_molar_mass(formula_input.upper())
                        if mass:
                            st.success(f"**Molar Mass:** {mass:.4f} g/mol")
                            with st.expander("📋 Composition Details"):
                                for elem in elements:
                                    st.write(elem)
                        else:
                            st.error(f"Error: {elements}")
            
            st.markdown("---")
            st.markdown("**Examples:**")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.code("H2O → 18.015 g/mol")
            with col2:
                st.code("C6H12O6 → 180.156 g/mol")
            with col3:
                st.code("NaCl → 58.44 g/mol")
        
        # Equation Balancer
        with tool_tabs[1]:
            st.markdown("##### Balance Chemical Equations")
            equation_input = st.text_input("Chemical Equation", placeholder="e.g., H2 + O2 -> H2O", key="balance_eq")
            
            if st.button("Balance", key="balance_btn"):
                if equation_input:
                    with st.spinner("Balancing..."):
                        balanced = balance_equation(equation_input)
                        if balanced:
                            st.success(f"**Balanced Equation:**\n\n{balanced}")
                        else:
                            st.info("Try one of these examples:")
                            st.code("""
                            H2 + O2 -> H2O
                            CH4 + O2 -> CO2 + H2O
                            Fe + O2 -> Fe2O3
                            C6H12O6 + O2 -> CO2 + H2O
                            """)
        
        # PubChem Search
        with tool_tabs[2]:
            st.markdown("##### Search PubChem Database")
            compound_search = st.text_input("Compound Name", placeholder="e.g., aspirin, caffeine, glucose", key="pubchem_search")
            
            if st.button("Search", key="search_pubchem_btn"):
                if compound_search:
                    with st.spinner("Searching PubChem..."):
                        result = search_pubchem(compound_search)
                        if result:
                            st.success("✅ Results found!")
                            col1, col2 = st.columns(2)
                            with col1:
                                st.metric("Formula", result['formula'])
                                st.metric("InChIKey", result['inchikey'][:14] + "...")
                            with col2:
                                st.metric("Molecular Weight", f"{result['weight']} g/mol")
                            st.markdown(f"[View on PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/{compound_search})")
                        else:
                            st.warning("Compound not found. Try a different name or check spelling.")
    
    # Tab 4: Chemistry Chat
    with tab4:
        st.markdown("#### 💬 Chemistry Assistant")
        
        # Initialize chat history
        if "messages" not in st.session_state:
            st.session_state.messages = [
                {"role": "assistant", "content": "👋 Hello! I'm Ojo Chem Pro. Ask me about molecules, symmetry, equations, or chemistry concepts!"}
            ]
        
        # Display chat history
        chat_container = st.container()
        with chat_container:
            for message in st.session_state.messages:
                if message["role"] == "user":
                    st.markdown(f'<div class="chat-message user-message" style="float: right; margin-left: auto;">{message["content"]}</div>', unsafe_allow_html=True)
                else:
                    st.markdown(f'<div class="chat-message bot-message" style="float: left;">{message["content"]}</div>', unsafe_allow_html=True)
                st.markdown("<br clear='both'>", unsafe_allow_html=True)
        
        # Chat input
        user_query = st.chat_input("Ask me about chemistry...")
        
        if user_query:
            # Add user message
            st.session_state.messages.append({"role": "user", "content": user_query})
            
            # Get response
            with st.spinner("Thinking..."):
                response = get_chatbot_response(user_query, safe_mode)
            
            # Add assistant response
            st.session_state.messages.append({"role": "assistant", "content": response})
            
            # Rerun to display new messages
            st.rerun()
        
        # Quick suggestions
        st.markdown("#### Quick Suggestions")
        cols = st.columns(4)
        suggestions = [
            "Show me benzene",
            "Symmetry of water",
            "Balance: H2 + O2 -> H2O",
            "Molar mass of glucose"
        ]
        
        for i, suggestion in enumerate(suggestions):
            with cols[i]:
                if st.button(suggestion, key=f"sugg_{i}"):
                    st.session_state.messages.append({"role": "user", "content": suggestion})
                    response = get_chatbot_response(suggestion, safe_mode)
                    st.session_state.messages.append({"role": "assistant", "content": response})
                    st.rerun()

if __name__ == "__main__":
    main()
