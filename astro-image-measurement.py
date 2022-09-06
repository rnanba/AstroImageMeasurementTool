#!/usr/bin/env python
import json
import math
import numpy as np
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle, Distance
from astropy.wcs import WCS
from astropy.wcs.utils import fit_wcs_from_points
from astroquery.simbad import Simbad
import svgwrite
from PIL import Image
import base64
from argparse import ArgumentParser

version = '0.1'

argparser = ArgumentParser(description='Calculate RA/DEC of target stars '\
                           'from the position of the star on the input image, '\
                           'and draw an annotated image.')
argparser.add_argument('ref_json', metavar='reference_stars.json',
                       help="data file of reference stars.")
argparser.add_argument('target_json', metavar='target_stars.json',
                       help="data file of target stars.")
argparser.add_argument('in_image', help="input image file (.jpg, .jpeg, .png)")
argparser.add_argument('out_image', metavar='out_image.svg',
                       help="output svg file of annotated image.")
argparser.add_argument('result_json', metavar='result.json', nargs='?',
                       help='output file (result data of target stars. '\
                       'if omitted, printed on stndard output.')
argparser.add_argument('--version', action='version',
                       version='%(prog)s ' + version)

args = argparser.parse_args()
ref_stars = None
obstime = None
target_stars = None

Simbad.add_votable_fields('pm')
Simbad.add_votable_fields('parallax')
Simbad.add_votable_fields('rv_value')

with open(args.ref_json, 'r', encoding='utf-8') as f:
    ref = json.load(f)
    obstime = Time(ref['obstime'])
    ref_stars = ref['stars']

with open(args.target_json, 'r', encoding='utf-8') as f:
    target_stars = json.load(f)

def ra_dec_to_angle(str_or_deg):
    if (type(str_or_deg) is str):
        return Angle(str_or_deg)
    else:
        return Angle(str_or_deg, u.deg)

def get_skycoord(name, obstime):
    rec = Simbad.query_object(name)[0]
    c =  SkyCoord(rec['RA'], rec['DEC'], unit=(u.hourangle, u.deg), frame="icrs",
                  pm_ra_cosdec=rec['PMRA']*u.mas/u.yr,
                  pm_dec=rec['PMDEC']*u.mas/u.yr,
                  distance=Distance(parallax=rec['PLX_VALUE']*u.mas),
                  radial_velocity=rec['RV_VALUE']*u.km/u.s,
                  obstime=Time('2000-01-01 11:58:55.816'))
    return c.apply_space_motion(new_obstime=Time(obstime))
    
ref_ra = []
ref_dec = []
ref_x = []
ref_y = []
proj_point = 'center'
for star in ref_stars:
    c = get_skycoord(star['name'], obstime)
    star['ra'] = c.ra.to_string(unit=u.hour)
    star['dec'] = c.dec.to_string(unit=u.degree)
    ref_ra.append(c.ra.deg)
    ref_dec.append(c.dec.deg)
    ref_x.append(star['x'])
    ref_y.append(star['y'])
    if ('proj_point' in star) and star['proj_point']:
        proj_point = c

target_xy = []
xy_index = []
target_radec = []
radec_index = []
for i, star in enumerate(target_stars):
    if 'x' in star:
        xy_index.append(i)
        target_xy.append([star['x'], star['y']])
    else:
        c = get_skycoord(star['name'], obstime)
        radec_index.append(i)
        star['ra'] = c.ra.to_string(unit=u.hour)
        star['dec'] = c.dec.to_string(unit=u.degree)
        target_radec.append((c.ra.deg, c.dec.deg))

stars = SkyCoord(ra=ref_ra, dec=ref_dec, frame="icrs", unit=u.deg)
pixels_x = np.array(ref_x)
pixels_y = np.array(ref_y)
wcs = fit_wcs_from_points(xy=[pixels_x, pixels_y],
                          world_coords=stars,
                          proj_point=proj_point)

radec_result = wcs.wcs_pix2world(target_xy, 0)
# for xy in target_xy:
#    print(wcs.pixel_to_world(xy[0], xy[1]))

for i, star_index in enumerate(xy_index):
    star = target_stars[star_index]
    star['ra'] = Angle(radec_result[i][0], u.deg).to_string(unit=u.hour)
    star['dec'] = Angle(radec_result[i][1], u.deg).to_string(unit=u.degree)

xy_result = wcs.wcs_world2pix(target_radec, 0)
for i, star_index in enumerate(radec_index):
    star = target_stars[star_index]
    star['x'] = xy_result[i][0]
    star['y'] = xy_result[i][1]

output = json.dumps(target_stars, indent=2, ensure_ascii=False)
if args.result_json:
    with open(args.result_json, 'w', encoding='utf-8') as out:
        out.write(output)
else:
    print(output)

image = Image.open(args.in_image)
w = image.width
h = image.height
corner_points = wcs.wcs_pix2world([[0, 0], [w, 0], [0, h], [w, h]], 0)

ra_min = Angle(min(list(map(lambda p: p[0], corner_points))), u.deg)
ra_max = Angle(max(list(map(lambda p: p[0], corner_points))), u.deg)
dec_min = Angle(min(list(map(lambda p: p[1], corner_points))), u.deg)
dec_max = Angle(max(list(map(lambda p: p[1], corner_points))), u.deg)

# print("RA: " + ra_min.to_string(unit=u.hour)
#       + " - " + ra_max.to_string(unit=u.hour))
# print("DEC: " + dec_min.to_string(unit=u.deg)
#       + " - " + dec_max.to_string(unit=u.deg))

ra1 = ra_min.hms
ra2 = ra_max.hms
dec1 = dec_min.dms
dec2 = dec_max.dms

ra_h_start = int(ra1.h)
ra_h_last = int(ra2.h)
ra_m_start = int(ra1.m)
ra_m_last = int(ra2.m)
ra_s_start = int(ra1.s)
ra_s_last = int(ra2.s)
ra_scale = []
ra_scale_deg = []
for ra_h in range(ra_h_start, ra_h_last + 1):
    if ra_h == ra_h_last:
        ra_m_stop = ra_m_last + 1
    else:
        ra_m_stop = 60
    for ra_m in range(ra_m_start, ra_m_stop):
        if (ra_h == ra_h_last) and (ra_m == ra_m_last):
            ra_s_stop = ra_s_last + 1
        else:
            ra_s_stop = 60
        for ra_s in range(ra_s_start, ra_s_stop):
            ra = (ra_h, ra_m, ra_s)
            ra_scale.append(ra)
            ra_scale_deg.append(Angle(ra, unit="hourangle").deg)
        ra_s_start = 0
    ra_m_start = 0

dec_d_start = int(dec1.d)
dec_d_last = int(dec2.d)
dec_m_start = int(dec1.m)
dec_m_last = int(dec2.m)
dec_s_start = int(dec1.s)
dec_s_last = int(dec2.s)
dec_scale = []
dec_scale_deg = []
for dec_d in range(dec_d_start, dec_d_last + 1):
    if dec_d == dec_d_last:
        dec_m_stop = dec_m_last + 1
    else:
        dec_m_stop = 60 if (dec_d >= 0) else 1
    for dec_m in range(dec_m_start, dec_m_stop):
        if (dec_d == dec_d_last) and (dec_m == dec_m_last):
            dec_s_stop = dec_s_last + 1
        else:
            dec_s_stop = 60 if (dec_d >= 0) else 1
        for dec_s in list(range(dec_s_start, dec_s_stop)):
            dec = (dec_d, dec_m, dec_s)
            dec_scale.append(dec)
            dec_scale_deg.append(Angle(dec, unit=u.deg).deg)
        dec_s_start = 0 if (dec_d >= 0) else -59
    dec_m_start = 0 if (dec_d >= 0) else -59

# draw output svg
drw = svgwrite.Drawing(args.out_image, size=(w, h))

# draw image
with open(args.in_image, 'rb') as f:
    mime = ''
    if args.in_image.upper().endswith('.JPG'):
        mime = 'image/jpeg'
    elif args.in_image.upper().endswith('.JPEG'):
        mime = 'image/jpeg'
    elif args.in_image.upper().endswith('.PNG'):
        mime = 'image/png'
    image_data = 'data:' + mime + ';base64,'
    image_data += base64.standard_b64encode(f.read()).decode()
    drw.add(drw.image(image_data, id='original_image'))

# draw grid
ra_lines_h = []
ra_lines_m = []
ra_lines_s = []
for i,ra in enumerate(ra_scale):
    coords = []
    for dec_deg in dec_scale_deg:
        coords.append((ra_scale_deg[i], dec_deg))
    line = wcs.world_to_pixel(SkyCoord(coords, frame="icrs", unit=u.deg))
    ra_h, ra_m, ra_s = ra
    # print((ra_h, ra_m, ra_s))
    if (ra_m == 0) and (ra_s == 0):
        ra_lines_h.append(line)
    elif ra_s == 0:
        ra_lines_m.append(line)
    else:
        ra_lines_s.append(line)

dec_lines_d = []
dec_lines_m = []
dec_lines_s = []

for i,dec in enumerate(dec_scale):
    coords = []
    for ra_deg in ra_scale_deg:
        coords.append((ra_deg, dec_scale_deg[i]))
    line = wcs.world_to_pixel(SkyCoord(coords, frame="icrs", unit=u.deg))
    dec_d, dec_m, dec_s = dec
    # print((dec_d, dec_m, dec_s))
    if (dec_m == 0) and (dec_s == 0):
        dec_lines_d.append(line)
    elif dec_s == 0:
        dec_lines_m.append(line)
    else:
        dec_lines_s.append(line)

def draw_lines(name, drw, lines, style):
    g = drw.g(id=name)
    for line in lines:
        points = []
        for i, x in enumerate(line[0]):
            points.append((x, line[1][i]))
        svg_line = drw.polyline(points)
        svg_line.update({ 'style' : style })
        g.add(svg_line)
    return g

grid_g = drw.g(id='grid')
ra_g = drw.g(id='ra')
ra_g.add(draw_lines('ra_h', drw, ra_lines_h, "stroke: #fcc; stroke-width: 2;"))
ra_g.add(draw_lines('ra_m', drw, ra_lines_m, "stroke: #ffc; stroke-width: 2;"))
ra_g.add(draw_lines('ra_s', drw, ra_lines_s, "stroke: #ccc; stroke-width: 1;"))
grid_g.add(ra_g)
dec_g = drw.g(id='dec')
dec_g.add(draw_lines('dec_d' ,drw, dec_lines_d, "stroke: #fcc; stroke-width: 2;"))
dec_g.add(draw_lines('dec_m', drw, dec_lines_m, "stroke: #ffc; stroke-width: 2;"))
dec_g.add(draw_lines('dec_s', drw, dec_lines_s, "stroke: #ccc; stroke-width: 1;"))
grid_g.add(dec_g)
drw.add(grid_g)

# draw center mark
# center = SkyCoord.from_pixel(w/2, h/2, wcs)
# min = Angle((0, 1, 0), unit=u.deg)
# c0 = center.directional_offset_by(Angle(0, unit=u.deg), min).to_pixel(wcs)
# c90 = center.directional_offset_by(Angle(90, unit=u.deg), min).to_pixel(wcs)
# c180 = center.directional_offset_by(Angle(180, unit=u.deg), min).to_pixel(wcs)
# c270 = center.directional_offset_by(Angle(270, unit=u.deg), min).to_pixel(wcs)

def to_xy(pixel):
    x = (pixel[0]).item()
    y = (pixel[1]).item()
    return (x, y)

# l1 = drw.line(start=to_xy(c0), end=to_xy(c180), stroke='red', stroke_width='1')
# l2 = drw.line(start=to_xy(c90), end=to_xy(c270), stroke='red', stroke_width='1')
# drw.add(l1)
# drw.add(l2)

def draw_markers(prefix, drw, stars, color, font_size):
    i = 0
    for star in stars:
        i += 1
        name = prefix + '_star_' + str(i)
        gs = drw.g(id=name)
        x = star['x']
        y = star['y']
        l1 = drw.line(start=(x - 4, y), end=(x + 4, y),
                      stroke=color, stroke_width='1')
        l2 = drw.line(start=(x, y - 4), end=(x, y + 4),
                      stroke=color, stroke_width='1')
        gc = drw.g(id=(name + '_crosshair'))
        gc.add(l1)
        gc.add(l2)
        gs.add(gc)
        gl = drw.g(id=(name + '_label'))
        t1 = drw.text(star['name'], x=[x+24], y=[y-(font_size * 2)],
                      fill=color, font_size=font_size,
                      id=(name + '_label_name'))
        t2 = drw.text(star['ra'], x=[x+24], y=[y-font_size],
                      fill=color, font_size=font_size,
                      id=(name + '_label_ra'))
        t3 = drw.text(star['dec'], x=[x+24], y=[y],
                      fill=color, font_size=font_size,
                      id=(name + '_label_dec'))
        gl.add(t1)
        gl.add(t2)
        gl.add(t3)
        gs.add(gl)
        drw.add(gs)

draw_markers('ref', drw, ref_stars, "#8f8", 48)
draw_markers('target', drw, target_stars, "#ff8", 24)

drw.save(pretty=True)
