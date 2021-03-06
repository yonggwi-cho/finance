import requests
import datetime
import pandas as pd
import csv
 
def daily_price_historical(symbol, comparison_symbol, limit=1, aggregate=1, exchange='Coincheck', allData='true'):
    url = 'https://min-api.cryptocompare.com/data/histoday?fsym={}&tsym={}&limit={}&aggregate={}&allData={}'\
            .format(symbol.upper(), comparison_symbol.upper(), limit, aggregate, allData)
    if exchange:
        url += '&e={}'.format(exchange)
    page = requests.get(url)
    data = page.json()['Data']
    df = pd.DataFrame(data)
    df['timestamp'] = [datetime.datetime.fromtimestamp(d) for d in df.time]
    return df
 
df = daily_price_historical("BTC","JPY")
df.to_csv("conincheck.csv")
