import requests
from bs4 import BeautifulSoup
import json
import os

def scrape_pesticides():
    url = 'http://www.bcpcpesticidecompendium.org/class_insecticides.html'
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')

    pesticides = []
    for li in soup.find_all('li'):
        a_tags = li.find_all('a')
        for a in a_tags:
            if 'href' in a.attrs and a['href'].endswith('.html'):  # Checking if the link is pointing to a detailed page
                name = a.text.strip()
                link = f"http://www.bcpcpesticidecompendium.org/{a['href']}"
                pesticide_details = scrape_pesticide_details(link)
                pesticide_details['name'] = name
                pesticides.append(pesticide_details)
    return pesticides

def scrape_pesticide_details(url):
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')

    cas_number = None
    inchi_key = None

    for tr in soup.find_all('tr'):
        th = tr.find('th')
        td = tr.find('td')
        if th and td:
            if 'CAS Reg. No.' in th.text:
                cas_number = td.text.strip()
            elif 'InChIKey' in th.text:
                inchi_key = td.text.strip()

    return {
        'cas_number': cas_number,
        'inchi_key': inchi_key
    }

def save_data(data, filename):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)

if __name__ == "__main__":
    data = scrape_pesticides()
    save_data(data, "C:/Users/Chris/Pesticides_Database/data/pesticides_data.json")
    print(f"Scraped data for {len(data)} pesticides")
