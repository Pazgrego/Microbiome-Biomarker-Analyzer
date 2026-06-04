import pandas as pd

def main():
    # Path to your microbiome CSV file
    file_path = "docs/OTU_Table_P37 microbiome to work on.csv"
    
    print("--- Starting data inspection ---")
    try:
        # Load the file using pandas
        df = pd.read_csv(file_path)
        
        print("\n📊 First 5 rows of the dataset:")
        print("-" * 40)
        print(df.head())
        print("-" * 40)
        
        print("\n📋 Detected column names:")
        print(df.columns.tolist())
        
        print(f"\n🔢 Total rows: {len(df)}, Total columns: {len(df.columns)}")

    except FileNotFoundError:
        print(f"❌ Error: The file was not found at path: {file_path}")
    except Exception as e:
        print(f"❌ Unexpected error while reading the file: {e}")

if __name__ == "__main__":
    main()