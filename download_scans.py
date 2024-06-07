import air_download.air_download as air
import argparse

args = argparse.Namespace()
args.cred_path = 'C:/Users/Fayyaz/Documents/GitHub/air_download/air_login.txt'
args.URL       = 'https://air.radiology.ucsf.edu/api/'
args.acc       = '10022969454'
args.profile   = 1
args.output    = '10022969454.zip'

air.main(args)