import requests
import io
import csv

def WISP_data_API_call(request_string, output_file=None):
    """ parse the request string to the cuurent WISP data service and optionally save output to file.
    See https://github.com/monocle-h2020/Access_to_WISPstation_API/blob/main/WISPcloud_API_documentation.ipynb for help on building request strings
    examples:
    WISP_data_API_call('REQUEST=GetDocumentation')
    WISP_data_API_call('REQUEST=ShowRecent')
    WISP_data_API_call('REQUEST=GetData&TIME=2018-05-11T09:00,2018-05-11T19:00&INSTRUMENT=WISPstation001')
    """
    username = 'public_access'
    password = 'WISPstation'
    Data_service = 'https://wispcloud.waterinsight.nl/api/query?SERVICE=Data&VERSION=1.0'
    url_string = f'{Data_service}&{request_string}'
    print(url_string)
    if output_file:
        datastring = io.StringIO(requests.get(url_string, auth=(username, password)).content.decode('utf-8'))
        print(datastring.read(), file=open(output_file, 'w'))
    else:
        print('confirming that no filename was provided to write output into')
    datastring = io.StringIO(requests.get(url_string, auth=(username, password)).content.decode('utf-8'))
    data = list(csv.DictReader(filter(lambda row: row[0] != '#', datastring), delimiter='\t'))
    return data