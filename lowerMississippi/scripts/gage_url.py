import pandas as pd
import re

df = pd.read_csv("xs/xs_catalog.csv")

def clean_gage(g):
    """Extract 8-digit number from any gage string or bytes literal."""
    if pd.isna(g):
        return ""
    s = str(g)
    # decode byte-like text such as "b'   07374000'"
    s = s.encode("utf-8", errors="ignore").decode("utf-8", errors="ignore") if isinstance(g, bytes) else s
    # find an 8-digit sequence
    m = re.search(r"\d{8}", s)
    return m.group(0) if m else ""

df["gages_clean"] = df["gages"].apply(clean_gage)
df["usgs_url"] = df["gages_clean"].apply(
    lambda g: f"https://waterdata.usgs.gov/monitoring-location/{g}" if g else ""
)

df.to_csv("xs/xs_catalog_with_urls.csv", index=False)
print("✅ Added cleaned 'usgs_url' column. Saved to xs/xs_catalog_with_urls.csv")
print(df[["gages", "gages_clean", "usgs_url"]].head(10))
