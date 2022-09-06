#!/usr/bin/env python
import sys
import json
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from argparse import ArgumentParser

argparser = ArgumentParser(description='get RA/DEC of target position.')
argparser.add_argument('stars_json', metavar='stars.json',
                       help="stars data file.")
argparser.add_argument('star_name_a', help="name of star A.")
argparser.add_argument('star_name_b', help="name of star B.")
args = argparser.parse_args()

stars = None
with open(args.stars_json, 'r', encoding='utf-8') as f:
    stars = json.load(f)

star_a = None
star_b = None
for star in stars:
    if star['name'] == args.star_name_a:
        star_a = star
    elif star['name'] == args.star_name_b:
        star_b = star

coord_a = SkyCoord(Angle(star_a['ra']).deg, Angle(star_a['dec']).deg, unit=u.deg)
coord_b = SkyCoord(Angle(star_b['ra']).deg, Angle(star_b['dec']).deg, unit=u.deg)
for_ra = SkyCoord(Angle(star_b['ra']).deg, Angle(star_a['dec']).deg, unit=u.deg)
for_dec = SkyCoord(Angle(star_a['ra']).deg, Angle(star_b['dec']).deg, unit=u.deg)

print(star_a['name'] + ": " + star_a['ra'] + " / " + star_a['dec'])
print(star_b['name'] + ": " + star_b['ra'] + " / " + star_b['dec'])
print("separation: " + coord_b.separation(coord_a).to_string(unit=u.deg) +
      " (" + coord_a.separation(for_ra).to_string(unit=u.deg) + " / " +
      coord_a.separation(for_dec).to_string(unit=u.deg) + ")")
