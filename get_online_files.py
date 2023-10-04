from bs4 import BeautifulSoup
import requests


def find_files(url):
    soup = BeautifulSoup(requests.get(url).text, features="lxml")

    hrefs = []

    for a in soup.find_all('a'):
        hrefs.append(a['href'])

    return hrefs

def get_image_urls(dates_list):
    file_path_list = []
    for date1 in dates_list:
        date1 = date1.replace("-", "/")
        url = "https://hesperia.gsfc.nasa.gov/hessidata/metadata/qlook_image/{}/".format(date1)

        list_of_links = find_files(url)

        for link in list_of_links:
            if "fsimg" in str(link):
                str_link = str(link)

                file_path_list.append(
                    "https://hesperia.gsfc.nasa.gov/hessidata/metadata/qlook_image/{}/{}".format(date1, str(link)))
    return file_path_list