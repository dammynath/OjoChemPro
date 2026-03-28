"""
ChemAssistant Pro - Advanced 3D Chemistry Lab with Symmetry Analysis
Streamlit Deployment Version
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
    page_title="ChemAssistant Pro",
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
        display: inline-block;
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
    
    .emergency-banner {
        background: #dc3545;
        color: white;
        padding: 0.5rem;
        text-align: center;
        border-radius: 5px;
        margin-bottom: 1rem;
        cursor: pointer;
        font-size: 0.85rem;
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
    
    div[data-testid="stVerticalBlock"] > div {
        gap: 0.5rem;
    }
</style>
""", unsafe_allow_html=True)

# ========== CHEMICAL DATABASE ==========
ATOMIC_MASSES = {
    'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81, 'C': 12.011,
    'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.18, 'Na': 22.99, 'Mg': 24.305,
    'Al': 26.982, 'Si': 28.086, 'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'K': 39.098,
    'Ca': 40.078, 'Fe': 55.845, 'Cu': 63.546, 'Zn': 65.38, 'Ag': 107.87, 'I': 126.90
}

# Molecule symmetry database
MOLECULE_SYMMETRY = {
    "water": {
        "name": "Water",
        "formula": "H₂O",
        "point_group": "C₂v",
        "elements": [
            {"type": "C₂", "axis": [0, 1, 0], "description": "2-fold rotation axis"},
            {"type": "σᵥ", "plane": [1, 0, 0], "description": "Vertical mirror plane"}
        ],
        "description": "Bent molecule with C₂v symmetry. Has a C₂ axis through the oxygen atom and two mirror planes."
    },
    "benzene": {
        "name": "Benzene",
        "formula": "C₆H₆",
        "point_group": "D₆h",
        "elements": [
            {"type": "C₆", "axis": [0, 1, 0], "description": "6-fold principal axis"},
            {"type": "C₂", "axis": [1, 0, 0], "description": "2-fold axes"},
            {"type": "σₕ", "plane": [0, 1, 0], "description": "Horizontal mirror plane"},
            {"type": "i", "center": [0, 0, 0], "description": "Inversion center"}
        ],
        "description": "Aromatic ring with D₆h symmetry. High symmetry with multiple rotation axes and mirror planes."
    },
    "methane": {
        "name": "Methane",
        "formula": "CH₄",
        "point_group": "Td",
        "elements": [
            {"type": "C₃", "axis": [1, 1, 1], "description": "3-fold axes through C-H bonds"},
            {"type": "C₂", "axis": [1, 0, 0], "description": "2-fold axes through midpoints"}
        ],
        "description": "Perfect tetrahedral molecule with Td symmetry. Highest symmetry among simple molecules."
    },
    "ammonia": {
        "name": "Ammonia",
        "formula": "NH₃",
        "point_group": "C₃v",
        "elements": [
            {"type": "C₃", "axis": [0, 1, 0], "description": "3-fold rotation axis"}
        ],
        "description": "Trigonal pyramidal with C₃v symmetry. Has a C₃ axis through the nitrogen atom."
    },
    "carbon_dioxide": {
        "name": "Carbon Dioxide",
        "formula": "CO₂",
        "point_group": "D∞h",
        "elements": [
            {"type": "C∞", "axis": [0, 1, 0], "description": "Infinite rotation axis"},
            {"type": "σₕ", "plane": [0, 1, 0], "description": "Horizontal mirror plane"},
            {"type": "i", "center": [0, 0, 0], "description": "Inversion center"}
        ],
        "description": "Linear molecule with D∞h symmetry. Has infinite rotation axis and inversion center."
    }
}

# Molecule 3D structure data (SMILES format for py3Dmol)
MOLECULE_STRUCTURES = {
    "water": "O",
    "benzene": "c1ccccc1",
    "methane": "C",
    "ammonia": "N",
    "carbon_dioxide": "O=C=O",
    "ethanol": "CCO",
    "glucose": "C1C(C(C(C(O1)O)O)O)O",
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O"
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
    """Balance chemical equations (simplified version)"""
    balanced_equations = {
        "h2+o2->h2o": "2H₂ + O₂ → 2H₂O",
        "ch4+o2->co2+h2o": "CH₄ + 2O₂ → CO₂ + 2H₂O",
        "c6h12o6+o2->co2+h2o": "C₆H₁₂O₆ + 6O₂ → 6CO₂ + 6H₂O",
        "na+cl2->nacl": "2Na + Cl₂ → 2NaCl",
        "n2+h2->nh3": "N₂ + 3H₂ → 2NH₃",
        "fe+o2->fe2o3": "4Fe + 3O₂ → 2Fe₂O₃"
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

# ========== 3D MOLECULE VISUALIZATION ==========
def create_3d_molecule(smiles, molecule_name):
    """Create 3D molecule visualization using py3Dmol"""
    if smiles not in MOLECULE_STRUCTURES.values():
        smiles = MOLECULE_STRUCTURES.get(molecule_name.lower(), "C")
    
    viewer = py3Dmol.view(width=500, height=400)
    viewer.addModel(smiles, 'smi')
    viewer.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    viewer.setBackgroundColor('0x1a1a2e')
    viewer.zoomTo()
    viewer.spin(True)
    
    return viewer

# ========== SYMMETRY VISUALIZATION ==========
def create_symmetry_visualization(molecule_name):
    """Create symmetry element visualization using Plotly"""
    sym_data = MOLECULE_SYMMETRY.get(molecule_name.lower())
    if not sym_data:
        return None
    
    fig = go.Figure()
    
    # Add sphere for molecule center
    fig.add_trace(go.Scatter3d(
        x=[0], y=[0], z=[0],
        mode='markers',
        marker=dict(size=8, color='white', symbol='circle'),
        name='Molecule Center'
    ))
    
    # Add symmetry elements
    colors = {'C₂': 'red', 'C₃': 'orange', 'C₆': 'yellow', 'σᵥ': 'green', 'σₕ': 'lightgreen', 'i': 'gold'}
    
    for element in sym_data["elements"]:
        if "axis" in element:
            # Draw rotation axis
            axis = element["axis"]
            fig.add_trace(go.Scatter3d(
                x=[-axis[0]*2, axis[0]*2],
                y=[-axis[1]*2, axis[1]*2],
                z=[-axis[2]*2, axis[2]*2],
                mode='lines',
                line=dict(color=colors.get(element["type"], 'red'), width=4),
                name=f"{element['type']} Axis"
            ))
        elif "plane" in element:
            # Draw mirror plane (as semi-transparent surface)
            plane = element["plane"]
            u, v = np.mgrid[-2:2:20j, -2:2:20j]
            if plane[1] != 0:  # Horizontal plane
                x = u
                z = v
                y = np.zeros_like(u)
            else:
                x = np.zeros_like(u)
                y = u
                z = v
            
            fig.add_trace(go.Surface(
                x=x, y=y, z=z,
                opacity=0.3,
                colorscale=[[0, 'green'], [1, 'lightgreen']],
                showscale=False,
                name=f"{element['type']} Plane"
            ))
        elif "center" in element:
            # Draw inversion center
            fig.add_trace(go.Scatter3d(
                x=[0], y=[0], z=[0],
                mode='markers',
                marker=dict(size=12, color='gold', symbol='circle-open', 
                           line=dict(color='gold', width=2)),
                name='Inversion Center'
            ))
    
    # Update layout
    fig.update_layout(
        title=f"Symmetry Elements of {sym_data['name']} (Point Group: {sym_data['point_group']})",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z",
            bgcolor='#1a1a2e',
            xaxis=dict(color='white'),
            yaxis=dict(color='white'),
            zaxis=dict(color='white')
        ),
        paper_bgcolor='#1a1a2e',
        plot_bgcolor='#1a1a2e',
        font=dict(color='white'),
        height=450
    )
    
    return fig

# ========== CHATBOT RESPONSE ENGINE ==========
def get_chatbot_response(query, safe_mode=True):
    """Generate response for user queries"""
    query_lower = query.lower()
    
    # Safety check
    hazardous_keywords = ['explosive', 'cyanide', 'mustard gas', 'weapon', 'toxic gas']
    if safe_mode and any(keyword in query_lower for keyword in hazardous_keywords):
        return "⚠️ **Safety Notice**: This topic involves hazardous materials and is blocked in Safe Mode. Please consult with a teacher or guardian."
    
    # Molecule visualization
    for molecule in MOLECULE_SYMMETRY.keys():
        if molecule in query_lower:
            return f"🔬 **{MOLECULE_SYMMETRY[molecule]['name']}**\n\n" \
                   f"**Formula:** {MOLECULE_SYMMETRY[molecule]['formula']}\n" \
                   f"**Point Group:** {MOLECULE_SYMMETRY[molecule]['point_group']}\n\n" \
                   f"**Description:** {MOLECULE_SYMMETRY[molecule]['description']}\n\n" \
                   f"**Symmetry Elements:**\n" + \
                   "\n".join([f"• {elem['description']}" for elem in MOLECULE_SYMMETRY[molecule]['elements']]) + \
                   f"\n\n💡 **Tip:** Check the 3D Viewer tab to see the molecule structure, and the Symmetry Analysis tab for 3D visualization of symmetry elements!"
    
    # Equation balancing
    if "balance" in query_lower or "->" in query:
        balanced = balance_equation(query)
        if balanced:
            return f"✅ **Balanced Equation:**\n\n{balanced}\n\n📝 **Steps:**\n1. Count atoms on each side\n2. Add coefficients\n3. Verify balance"
        return "🤔 I can help balance equations! Try examples like:\n• H2 + O2 -> H2O\n• CH4 + O2 -> CO2 + H2O\n• C6H12O6 + O2 -> CO2 + H2O"
    
    # Molar mass calculation
    if "molar mass" in query_lower or "molecular weight" in query_lower:
        # Extract formula
        import re
        formula_match = re.search(r'([A-Z][a-z]?\d*)+', query)
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
        return "Please specify a compound to search (e.g., 'Search PubChem for aspirin')"
    
    # General chemistry help
    return f"🧪 **Chemistry Assistant**\n\nI can help you with:\n\n" \
           f"• **3D Molecule Viewer** - View molecules in 3D\n" \
           f"• **Symmetry Analysis** - Visualize symmetry elements\n" \
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
        <h1 style="color: white; margin: 0;">🧪 ChemAssistant Pro</h1>
        <p style="color: white; margin: 0; opacity: 0.9;">Advanced 3D Chemistry Lab with Symmetry Analysis</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Emergency banner (minimized)
    with st.expander("🚨 Emergency Information", expanded=False):
        st.warning("""
        **Poison Control:** 1-800-222-1222  
        **CHEMTREC:** 1-800-424-9300  
        
        **First Aid:** Flush with water for 15 minutes, seek immediate medical attention.
        """)
    
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
        cols = st.columns(3)
        molecules = ["Water", "Benzene", "Methane", "Ammonia", "CO₂", "Ethanol"]
        
        for i, mol in enumerate(molecules):
            with cols[i % 3]:
                if st.button(f"🔬 {mol}", key=f"quick_{mol}", use_container_width=True):
                    st.session_state['selected_molecule'] = mol.lower()
                    st.session_state['active_tab'] = "Viewer"
                    st.rerun()
        
        st.divider()
        
        # Mode info
        st.markdown("### ℹ️ About")
        st.info("""
        **ChemAssistant Pro** helps you explore:
        - 3D molecular structures
        - Symmetry elements and point groups
        - Chemical calculations
        - Equation balancing
        """)
        
        # Stats
        st.divider()
        st.markdown(f"**Molecules:** {len(MOLECULE_SYMMETRY)}")
        st.markdown(f"**Elements:** {len(ATOMIC_MASSES)}")
        st.markdown(f"**Last Updated:** {datetime.now().strftime('%Y-%m-%d')}")
    
    # Main content area with tabs
    tab1, tab2, tab3, tab4 = st.tabs(["🔬 3D Viewer", "🔄 Symmetry Analysis", "⚙️ Tools", "💬 Chemistry Chat"])
    
    # Tab 1: 3D Viewer
    with tab1:
        col1, col2 = st.columns([2, 1])
        
        with col1:
            st.markdown("#### 3D Molecular Structure")
            
            # Molecule selector
            molecule_name = st.selectbox(
                "Select Molecule",
                options=list(MOLECULE_STRUCTURES.keys()),
                format_func=lambda x: x.title(),
                key="molecule_selector"
            )
            
            if 'selected_molecule' in st.session_state:
                molecule_name = st.session_state['selected_molecule']
            
            # Display 3D molecule
            smiles = MOLECULE_STRUCTURES.get(molecule_name, "C")
            viewer = create_3d_molecule(smiles, molecule_name)
            showmol(viewer, height=450, width=600)
            
            # Molecule info
            if molecule_name in MOLECULE_SYMMETRY:
                sym = MOLECULE_SYMMETRY[molecule_name]
                st.info(f"""
                **{sym['name']}** ({sym['formula']})  
                **Point Group:** {sym['point_group']}  
                **Description:** {sym['description']}
                """)
        
        with col2:
            st.markdown("#### 🎮 Viewer Controls")
            st.markdown("""
            - **🖱️ Drag** - Rotate view
            - **🖱️ Right-click drag** - Pan
            - **🖱️ Scroll** - Zoom in/out
            - **🔄 Auto-rotate** - Animation enabled
            """)
            
            st.markdown("#### 📊 Molecule Info")
            st.markdown(f"""
            - **Formula:** {MOLECULE_SYMMETRY.get(molecule_name, {}).get('formula', 'N/A')}
            - **Atoms:** Varies by molecule
            - **Bonds:** Varies by molecule
            """)
            
            if molecule_name in MOLECULE_SYMMETRY:
                with st.expander("🔍 Symmetry Elements"):
                    for elem in MOLECULE_SYMMETRY[molecule_name]['elements']:
                        st.markdown(f"- **{elem['type']}**: {elem['description']}")
    
    # Tab 2: Symmetry Analysis
    with tab2:
        st.markdown("#### 🔄 Symmetry Elements Visualization")
        
        col1, col2 = st.columns([1.5, 1])
        
        with col1:
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
                st.info("Symmetry visualization loading...")
        
        with col2:
            sym_data = MOLECULE_SYMMETRY.get(sym_molecule, {})
            st.markdown(f"### {sym_data.get('name', 'N/A')}")
            st.markdown(f"**Formula:** {sym_data.get('formula', 'N/A')}")
            st.markdown(f"**Point Group:** `{sym_data.get('point_group', 'N/A')}`")
            
            st.markdown("#### Symmetry Elements")
            for elem in sym_data.get('elements', []):
                st.markdown(f"- **{elem['type']}**: {elem['description']}")
            
            st.markdown("#### Symmetry Guide")
            with st.expander("📖 Understanding Symmetry"):
                st.markdown("""
                **Rotation Axes (Cₙ):** Axis around which molecule rotates to appear identical.
                
                **Mirror Planes (σ):** Planes that reflect molecule onto itself.
                
                **Inversion Center (i):** Point where (x,y,z) → (-x,-y,-z).
                
                **Point Groups:** Classification of molecular symmetry.
                """)
    
    # Tab 3: Tools
    with tab3:
        st.markdown("#### ⚙️ Chemistry Tools")
        
        tool_tabs = st.tabs(["📊 Molar Mass", "⚖️ Equation Balancer", "🔍 PubChem Search"])
        
        # Molar Mass Calculator
        with tool_tabs[0]:
            st.markdown("##### Calculate Molar Mass")
            formula_input = st.text_input("Chemical Formula", placeholder="e.g., H2O, C6H12O6, NaCl")
            
            if st.button("Calculate", key="calc_molar"):
                if formula_input:
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
            st.code("""
            Water (H2O): 18.015 g/mol
            Glucose (C6H12O6): 180.156 g/mol
            Sodium Chloride (NaCl): 58.44 g/mol
            """)
        
        # Equation Balancer
        with tool_tabs[1]:
            st.markdown("##### Balance Chemical Equations")
            equation_input = st.text_input("Chemical Equation", placeholder="e.g., H2 + O2 -> H2O")
            
            if st.button("Balance", key="balance_eq"):
                if equation_input:
                    balanced = balance_equation(equation_input)
                    if balanced:
                        st.success(f"**Balanced Equation:**\n\n{balanced}")
                    else:
                        st.info("Try one of these examples:")
                        st.code("""
                        H2 + O2 -> H2O
                        CH4 + O2 -> CO2 + H2O
                        Fe + O2 -> Fe2O3
                        """)
        
        # PubChem Search
        with tool_tabs[2]:
            st.markdown("##### Search PubChem Database")
            compound_search = st.text_input("Compound Name", placeholder="e.g., aspirin, caffeine, glucose")
            
            if st.button("Search", key="search_pubchem"):
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
                {"role": "assistant", "content": "👋 Hello! I'm ChemAssistant Pro. Ask me about molecules, symmetry, equations, or chemistry concepts!"}
            ]
        
        # Display chat history
        for message in st.session_state.messages:
            if message["role"] == "user":
                st.markdown(f'<div class="chat-message user-message" style="float: right;">{message["content"]}</div>', unsafe_allow_html=True)
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