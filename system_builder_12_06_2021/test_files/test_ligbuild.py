def collect_smi_data(input_smi_file):
    """
    Reads in a .smi file with format: SMILES_STRING LIGAND_NAME
    
    Args:
        input_smi_file (str): Path to the .smi file
        
    SMI file format:
        c1ccc(c(c1)[C@@H]2CCN(C2)c3c4cc[nH]c4ncn3)F CACHE3HI_1715_9_1
        c1ccc(c(c1)[C@H]2CCN(C2)c3c4cc[nH]c4ncn3)F CACHE3HI_1715_9_2
    """
    
    import os
    
    if not os.path.exists(input_smi_file):
        raise FileNotFoundError(f"SMI file not found: {input_smi_file}")
    
    smiles_list = []
    names_list = []
    
    print(f"Reading SMI file: {input_smi_file}")
    
    with open(input_smi_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            
            # Split by whitespace (space or tab)
            parts = line.split()
            
            if len(parts) < 2:
                print(f"Warning: Line {line_num} doesn't have both SMILES and name. Skipping: {line[:50]}")
                continue
            
            # First part is SMILES, everything else is the name
            smiles = parts[0]
            name = ' '.join(parts[1:])
            
            smiles_list.append(smiles)
            names_list.append(name)
    
    print(f"Successfully loaded {len(smiles_list)} ligands from {input_smi_file}")
    
    if len(smiles_list) == 0:
        raise ValueError(f"No valid ligands found in {input_smi_file}")
    
    return smiles_list, names_list


# ============================================================================
# TEST SCRIPT - Run this to verify it works!
# ============================================================================

if __name__ == "__main__":
    print("Testing collect_smi_data()...")
    print("=" * 60)
    
    # Create a test .smi file
    test_content = """# Test SMI file
c1ccc(c(c1)[C@@H]2CCN(C2)c3c4cc[nH]c4ncn3)F CACHE3HI_1715_9_1
c1ccc(c(c1)[C@H]2CCN(C2)c3c4cc[nH]c4ncn3)F CACHE3HI_1715_9_2

# Empty line above is ignored
CN1CCN(CC1)c2ccc(cc2)NC(=O)c3ccc(cc3)C(C)(C)C test_ligand_3
"""
    
    with open('test.smi', 'w') as f:
        f.write(test_content)
    
    # Test the function
    try:
        smiles, names = collect_smi_data('test.smi')
        
        print(f"\nSuccessfully loaded {len(smiles)} ligands")
        print("\nLigands found:")
        for s, n in zip(smiles, names):
            print(f"  {n}")
            print(f"    SMILES: {s[:60]}...")
            print()
        
        print("=" * 60)
        print("SUCCESS! The function works correctly.")
        print("\nNext steps:")
        print("1. Copy this function into your ligand_builder_v0_3.py")
        print("2. Test with your actual .smi file")
        print("3. Commit: git commit -m 'feat: add SMI file support'")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
    
    finally:
        # Cleanup
        import os
        if os.path.exists('test.smi'):
            os.remove('test.smi')
