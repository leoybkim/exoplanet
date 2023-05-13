import requests
import pandas as pd
from tqdm import tqdm

# Set the API URL and parameters
url = "https://exoplanets.nasa.gov/api/v1/planets"
params = {"order": "display_name+asc", "per_page": "50", "page": "0", "search": ""}

# Initialize an empty list to store the data
data = []

# Loop through all pages and append the data to the list
for page in tqdm(range(107)):
    params["page"] = page
    response = requests.get(url, params=params)
    if response.ok:
        response_data = response.json()
        data += response_data["items"]
    else:
        print("Error: Could not retrieve data from API.")

# Convert the data list to a dataframe
df = pd.DataFrame(data)

# Save the dataframe to a csv file
df.to_csv("nasa_exoplanet_catalog.csv", index=False)
