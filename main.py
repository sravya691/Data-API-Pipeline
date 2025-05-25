import streamlit as st
from SPARQLWrapper import SPARQLWrapper, JSON
from rapidfuzz.fuzz import ratio
import urllib.parse
import re

# Must be the very first Streamlit call
st.set_page_config(page_title="MetaNetX Compound Search", layout="wide")

# Query parameters - FIXED
query_params = st.query_params
compound_param = query_params.get("compound")
view_param = query_params.get("view")

# Helper Functions
def extract_compounds_from_expression(expression):
    """Extract individual compound names from chemical expressions, equations, or pathways"""
    # Common separators in chemical expressions
    separators = [
        '→', '⟶', '->', '⇒', '⇆', '⇌', '↔',  # arrows
        '+', '＋',  # plus signs
        '⟷', '↔',  # equilibrium arrows
        '|', '/', '\\',  # alternative separators
        ',', ';',  # list separators
        '⟨', '⟩', '[', ']', '(', ')',  # brackets (split around them)
    ]
    
    # Replace all separators with a common delimiter
    cleaned_expression = expression
    for sep in separators:
        cleaned_expression = cleaned_expression.replace(sep, '|')
    
    # Split and clean compounds
    compounds = []
    for part in cleaned_expression.split('|'):
        part = part.strip()
        if part and len(part) > 0:
            # Remove stoichiometric coefficients (numbers at the beginning)
            part = re.sub(r'^\d+\s*', '', part)
            # Remove common prefixes/suffixes that aren't part of compound names
            part = re.sub(r'\s*(aq|s|l|g)\s*$', '', part)  # phase indicators
            part = part.strip()
            if part and len(part) > 1:  # Avoid single characters
                compounds.append(part)
    
    # Remove duplicates while preserving order
    seen = set()
    unique_compounds = []
    for comp in compounds:
        if comp.lower() not in seen:
            seen.add(comp.lower())
            unique_compounds.append(comp)
    
    return unique_compounds

def normalize_compound_name(name):
    """Normalize compound names for better database matching"""
    if not name:
        return ""
    
    # Handle Unicode subscripts and superscripts
    subscript_map = {'₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4', 
                     '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9'}
    superscript_map = {'⁰': '0', '¹': '1', '²': '2', '³': '3', '⁴': '4', 
                       '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9',
                       '⁺': '+', '⁻': '-'}
    
    normalized = name
    for unicode_char, ascii_char in {**subscript_map, **superscript_map}.items():
        normalized = normalized.replace(unicode_char, ascii_char)
    
    # Additional normalizations
    normalized = normalized.replace('–', '-').replace('—', '-')  # different dashes
    normalized = normalized.replace(''', "'").replace(''', "'")  # different quotes
    
    return normalized.strip()

def is_likely_equation(text):
    """Check if the input looks like a chemical equation rather than a single compound"""
    equation_indicators = ['→', '⟶', '->', '⇒', '⇆', '⇌', '↔', '+', '⟷']
    return any(indicator in text for indicator in equation_indicators)

def get_pubchem_img_url(name):
    try:
        encoded_name = urllib.parse.quote(normalize_compound_name(name))
        return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/PNG"
    except Exception:
        return None

def get_external_refs(compound_name):
    normalized_name = normalize_compound_name(compound_name)
    query = f"""
    PREFIX mnx: <https://rdf.metanetx.org/schema/>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

    SELECT ?metabolite ?xref
    WHERE {{
        ?metabolite a mnx:CHEM .
        ?metabolite rdfs:comment ?comment .
        FILTER(LCASE(?comment) = "{normalized_name.lower()}")
        ?metabolite mnx:chemXref ?xref
    }}
    """
    sparql = SPARQLWrapper("https://rdf.metanetx.org/sparql")
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    try:
        results = sparql.query().convert()
        return [r["xref"]["value"] for r in results["results"]["bindings"]]
    except Exception as e:
        st.error(f"External references query failed: {e}")
        return []

def make_clickable_xref(xref):
    if ":" not in xref:
        return xref
    db, identifier = xref.split(":", 1)
    links = {
        "CHEBI": f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{identifier}",
        "KEGG": f"https://www.kegg.jp/dbget-bin/www_bget?cpd:{identifier}",
        "HMDB": f"https://hmdb.ca/metabolites/{identifier}",
        "PubChem": f"https://pubchem.ncbi.nlm.nih.gov/compound/{identifier}",
        "MetaCyc": f"https://metacyc.org/compound?orgid=META&id={identifier}"
    }
    url = links.get(db)
    return f"[{xref}]({url})" if url else xref

def display_external_refs(compound_name):
    xrefs = get_external_refs(compound_name)
    if xrefs:
        for xref in xrefs:
            st.markdown(f"- {make_clickable_xref(xref)}")
    else:
        st.markdown("No references found.")

# Special reference-only page
if compound_param and view_param == "refs":
    st.title(f"External References for: {compound_param}")
    display_external_refs(compound_param)
    st.markdown("[Back to search](/)", unsafe_allow_html=True)
    st.stop()

def run_query(compound_name, use_contains=False):
    endpoint_url = "https://rdf.metanetx.org/sparql"
    sparql = SPARQLWrapper(endpoint_url)
    normalized_name = normalize_compound_name(compound_name)
    
    # Escape quotes in the compound name for SPARQL
    escaped_name = normalized_name.replace('"', '\\"').replace("'", "\\'")
    
    filter_clause = (
        f'FILTER(CONTAINS(LCASE(?comment), "{escaped_name.lower()}"))'
        if use_contains else
        f'FILTER(LCASE(?comment) = "{escaped_name.lower()}")'
    )
    query = f"""
    PREFIX mnx: <https://rdf.metanetx.org/schema/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

    SELECT ?metabolite ?label ?comment ?reference ?formula ?charge ?inchi ?inchikey ?smiles
    WHERE {{
        ?metabolite a mnx:CHEM .
        ?metabolite rdfs:label ?label .
        ?metabolite rdfs:comment ?comment .
        {filter_clause}
        ?metabolite mnx:chemRefer ?reference .
        OPTIONAL {{ ?metabolite mnx:formula  ?formula }}
        OPTIONAL {{ ?metabolite mnx:charge   ?charge }}
        OPTIONAL {{ ?metabolite mnx:inchi    ?inchi }}
        OPTIONAL {{ ?metabolite mnx:inchikey ?inchikey }}
        OPTIONAL {{ ?metabolite mnx:smiles   ?smiles }}
    }}
    """
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    try:
        results = sparql.query().convert()
        return results.get("results", {}).get("bindings", [])
    except Exception as e:
        st.error(f"SPARQL query failed: {e}")
        return []

def fetch_data(compound_name, match_type="Both", partial_limit=5):
    results = {"exact": [], "partial": []}
    if match_type in ["Exact", "Both"]:
        results["exact"] = run_query(compound_name, use_contains=False)
    if match_type in ["Partial", "Both"]:
        partials = run_query(compound_name, use_contains=True)
        normalized_input = normalize_compound_name(compound_name).lower()
        partials = [
            r for r in partials
            if normalize_compound_name(r.get("comment", {}).get("value", "")).lower() != normalized_input
        ]
        scored_partials = []
        for r in partials:
            name = r.get("comment", {}).get("value", "")
            sim = ratio(normalize_compound_name(name).lower(), normalized_input)
            scored_partials.append((sim, r))
        scored_partials.sort(reverse=True, key=lambda x: x[0])
        results["partial"] = scored_partials[:partial_limit]
    return results

def display_compound_results(results, compound_name, is_multi_search=False):
    """Display results for a single compound"""
    if results["exact"]:
        header = f"Exact Match for '{compound_name}'" if is_multi_search else "Exact Match"
        st.subheader(header)
        display_results_list(results["exact"], show_similarity=False)
    elif is_multi_search or results.get("searched_exact", False):
        if is_multi_search:
            st.info(f"No exact match found for '{compound_name}'")

    if results["partial"]:
        header = f"Top Partial Matches for '{compound_name}'" if is_multi_search else "Top Partial Matches"
        st.subheader(header)
        display_results_list(results["partial"], show_similarity=True)

def display_results_list(results_list, show_similarity=False):
    """Display a list of compound results"""
    for i, item in enumerate(results_list, 1):
        r = item[1] if show_similarity else item
        score = item[0] if show_similarity else None
        
        st.markdown("---")
        comment = r.get('comment', {}).get('value', 'N/A')
        
        try:
            encoded_comment = urllib.parse.quote(comment, safe='')
        except:
            encoded_comment = comment
            
        label = r.get('label', {}).get('value', 'N/A')
        uri = r.get('metabolite', {}).get('value', '')
        ref = r.get('reference', {}).get('value', '')
        formula = r.get('formula', {}).get('value', 'N/A')
        charge = r.get('charge', {}).get('value', 'N/A')
        inchi = r.get('inchi', {}).get('value', 'N/A')
        inchikey = r.get('inchikey', {}).get('value', 'N/A')
        smiles = r.get('smiles', {}).get('value', 'N/A')
        
        col1, col2 = st.columns([2, 1])
        with col1:
            if show_similarity:
                st.markdown(f"**#{i} — Similarity: {score:.1f}%**")
            st.markdown(f"**Label:** `{comment}`")
            st.markdown(f"**MNX ID:** {label}")
            st.markdown(f"**URI:** [{uri}]({uri})")
            st.markdown(f"**Reference:** [{ref}]({ref})")
            st.markdown(f"**Formula:** `{formula}`")
            st.markdown(f"**Charge:** `{charge}`")
            st.markdown(f"**InChI:** `{inchi}`")
            st.markdown(f"**InChIKey:** `{inchikey}`")
            st.markdown(f"**SMILES:** `{smiles}`")
            st.markdown(f"**External References:** [View references](?compound={encoded_comment}&view=refs)")
        with col2:
            img_url = get_pubchem_img_url(comment)
            if img_url:
                try:
                    st.image(img_url, width=150)
                    st.markdown(f"[Image source: PubChem]({img_url})")
                except:
                    st.markdown("Image not available")

# Interface
st.title("MetaNetX Compound Explorer")

# Input handling
compound_input = st.text_input("Enter compound name or chemical equation:")

match_type = st.selectbox("Select match type:", ["Exact", "Partial", "Both"])
partial_limit = st.number_input("Limit partial matches:", min_value=1, max_value=50, value=5) if match_type in ["Partial", "Both"] else 5

if st.button("Search"):
    if not compound_input.strip():
        st.warning("Please enter a compound name or chemical expression.")
    else:
        # Check if input looks like an equation/expression
        if is_likely_equation(compound_input):
            st.info("Detected chemical equation/expression. Searching for individual compounds...")
            compounds = extract_compounds_from_expression(compound_input)
            
            if compounds:
                st.markdown(f"**Extracted compounds:** {', '.join(compounds)}")
                
                # Search each compound
                for compound in compounds:
                    if len(compound) > 1:  # Skip very short strings
                        st.markdown("---")
                        results = fetch_data(compound, match_type, partial_limit)
                        results["searched_exact"] = match_type in ["Exact", "Both"]
                        display_compound_results(results, compound, is_multi_search=True)
            else:
                st.warning("Could not extract valid compound names from the expression.")
        else:
            # Single compound search
            results = fetch_data(compound_input, match_type, partial_limit)
            display_compound_results(results, compound_input, is_multi_search=False)
            
            if not results["exact"] and not results["partial"]:
                st.info("No matches found. Try using partial matching or check your compound name.")

st.markdown("---")
st.caption("Data from MetaNetX | External database links included for cross-referencing")
