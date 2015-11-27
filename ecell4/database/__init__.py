# from ecell4.core import Species
import urllib.request

def description(obj, database="uniprot"):
    if database == "uniprot":
        return description_uniprot(obj)
    return None

def description_uniprot(obj):
    # if isinstance(obj, Species) and obj.has_attribute("uniprot.id"):
    #     uid = obj.get_attribute("uniprot.id")
    # elif isinstance(obj, str):
    if isinstance(obj, str):
        uid = obj
    else:
        return None

    url = 'http://www.uniprot.org/uniprot/{}.txt'.format(uid)
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as response:
        data = response.read()
    return data.decode('utf-8')


if __name__ == "__main__":
    # sp = Species("MinD")
    # sp.set_attribute("uniprot.id", "P0AEZ3")
    # print(description(sp, database="uniprot")) 
    print(description("P0AEZ3", database="uniprot"))