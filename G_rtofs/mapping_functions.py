# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 12:10:58 2018

A series of functions for mapping between spheres and planes. Useful for working with
geophysical data.
@author: pnola
"""

import numpy

#Cotangent
def cot(th):
    return 1.0/numpy.tan(th)

#Secant
def sec(th):
    return 1.0/numpy.cos(th)

#Degrees to radians
def deg2rad(deg):
    return deg*numpy.pi/180.0

#Radians to degress
def rad2deg(rad):
    return rad*180.0/numpy.pi

#This function takes a spherical grid and maps it to a plane using a
# Lambert conformal projection.    
def lonlat2km(ref_lon,ref_lat,lon,lat,std_lat1=30,std_lat2=60):
    stdlat1  =  deg2rad(std_lat1)
    stdlat2  =  deg2rad(std_lat2)
    R=6371
    ref_lon = deg2rad(ref_lon)
    ref_lat = deg2rad(ref_lat)
    lon = deg2rad(lon)
    lat = deg2rad(lat)
    n = numpy.log(numpy.cos(stdlat1)*sec(stdlat2)) / numpy.log(numpy.tan(0.25*numpy.pi+0.5*stdlat2)*cot(0.25*numpy.pi+0.5*stdlat1))
    F=(numpy.cos(stdlat1)*(numpy.tan(0.25*numpy.pi+0.5*stdlat1)**n))/n
    p0 = R*F*(cot(0.25*numpy.pi+0.5*ref_lat)**n)
    p = R*F*(cot(0.25*numpy.pi+0.5*lat)**n)
    th = n*(lon-ref_lon)
    x=p*numpy.sin(th)
    y=p0-p*numpy.cos(th)
    return x,y

#This function takes a spherical grid and maps it to a plane using a
# Lambert conformal projection.    
def lonlat2m(ref_lon,ref_lat,lon,lat,std_lat1=30,std_lat2=60):
    stdlat1  =  deg2rad(std_lat1)
    stdlat2  =  deg2rad(std_lat2)
    R=6371
    ref_lon = deg2rad(ref_lon)
    ref_lat = deg2rad(ref_lat)
    lon = deg2rad(lon)
    lat = deg2rad(lat)
    n = numpy.log(numpy.cos(stdlat1)*sec(stdlat2)) / numpy.log(numpy.tan(0.25*numpy.pi+0.5*stdlat2)*cot(0.25*numpy.pi+0.5*stdlat1))
    F=(numpy.cos(stdlat1)*(numpy.tan(0.25*numpy.pi+0.5*stdlat1)**n))/n
    p0 = R*F*(cot(0.25*numpy.pi+0.5*ref_lat)**n)
    p = R*F*(cot(0.25*numpy.pi+0.5*lat)**n)
    th = n*(lon-ref_lon)
    x=p*numpy.sin(th)
    y=p0-p*numpy.cos(th)
    return x*1000,y*1000

#This function takes a cartesian grid and maps it to a sphere using an inverse
# Lambert conformal projection.
def km2lonlat(ref_lon,ref_lat,x,y,std_lat1=30,std_lat2=60):
    stdlat1  =  deg2rad(std_lat1)
    stdlat2  =  deg2rad(std_lat2)
    R=6371
    ref_lon = deg2rad(ref_lon)
    ref_lat = deg2rad(ref_lat)
    n = numpy.log(numpy.cos(stdlat1)*sec(stdlat2)) / numpy.log(numpy.tan(0.25*numpy.pi+0.5*stdlat2)*cot(0.25*numpy.pi+0.5*stdlat1))
    F=(numpy.cos(stdlat1)*(numpy.tan(0.25*numpy.pi+0.5*stdlat1)**n))/n
    p0 = R*F*(cot(0.25*numpy.pi+0.5*ref_lat)**n)
    p = numpy.sign(n)*numpy.sqrt(x**2+(p0-y)**2)
    th = numpy.arctan(x/(p0-y))
    lon = th/n + ref_lon
    lat = 2*numpy.arctan((R*F/p)**(1/n))-numpy.pi/2
    lon = rad2deg(lon)
    lat = rad2deg(lat)
    return lon,lat

def m2lonlat(ref_lon,ref_lat,x,y,std_lat1=30,std_lat2=60):
    x = x/1000
    y = y/1000
    stdlat1  =  deg2rad(std_lat1)
    stdlat2  =  deg2rad(std_lat2)
    R=6371
    ref_lon = deg2rad(ref_lon)
    ref_lat = deg2rad(ref_lat)
    n = numpy.log(numpy.cos(stdlat1)*sec(stdlat2)) / numpy.log(numpy.tan(0.25*numpy.pi+0.5*stdlat2)*cot(0.25*numpy.pi+0.5*stdlat1))
    F=(numpy.cos(stdlat1)*(numpy.tan(0.25*numpy.pi+0.5*stdlat1)**n))/n
    p0 = R*F*(cot(0.25*numpy.pi+0.5*ref_lat)**n)
    p = numpy.sign(n)*numpy.sqrt(x**2+(p0-y)**2)
    th = numpy.arctan(x/(p0-y))
    lon = th/n + ref_lon
    lat = 2*numpy.arctan((R*F/p)**(1/n))-numpy.pi/2
    lon = rad2deg(lon)
    lat = rad2deg(lat)
    return lon,lat