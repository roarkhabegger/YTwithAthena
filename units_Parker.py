#!/usr/bin/python
# Class for unit systems
# RSH 2/8/21
import astropy.units as u
import astropy.constants as c
import numpy as np

class Parker:
    def __init__(self,alpha=0.0,vx=0,gravH=100,scaleH=100,d=10**(-24),p=10**(-12)):
        a0              = 10**6 *u.cm/u.s
        g0              = 4*10**(-9) * u.cm /u.s**2
        rho0            = 10**(-24) *u.g/u.cm**3

        self.vs              = (a0).to('cm/s')
        self.ls              = (a0**2*(1+alpha)/g0).to('pc')
        self.ms              = (rho0*self.ls**3).to('g')
        self.ts              = (self.ls/self.vs).to('Myr')
        self.es              = (rho0*self.vs**2).to('erg/cm^3')
        self.rhos            = (rho0).to('g/cm^3')
        self.ns              = (1.0/self.ls**3).to('1/cm^3')

        xVel = vx*u.cm/u.s
        Hg = self.ls #gravH*u.pc
        H  = self.ls #scaleH*u.pc
        rho = d*u.g/u.cm**3
        grav = g0
        pres = p*u.erg/u.cm**3

        self.vc              = (xVel/ self.vs).to('')
        self.H0c             = (H/ self.ls).to('')
        self.Hgc             = (Hg/ self.ls).to('')
        self.rhoc            = (rho /self.rhos).to('')
        self.pc              = (pres/self.es).to('')
        self.gc              = (grav/(self.ls/self.ts**2)).to('')


    def printComp(self):
        print("Computational Mas Density =",format(self.rhoc,'1.3e'))
        print("Computational Gas Press   =",format(self.pc,'1.3e'))
        print("Computational Grav Accel  =",format(self.gc,'1.3e'))
        print("Computational Scale Hgt   =",format(self.H0c,'1.3e'))
        print("Computational Grav Hgt    =",format(self.Hgc,'1.3e'))
        print("Computational xVelocity   =",format(self.vc,'1.3e'))
        return
    def printScale(self):
        print("Distance Scaling = ",format(self.ls,'1.3e'))
        print("Density Scaling  = ",format(self.rhos,'1.3e'))
        print("Time Scaling     = ",format(self.ts,'1.3e'))
        print("Velocity Scaling = ",format(self.vs,'1.3e'))
        print("Pressure Scaling = ",format(self.es,'1.3e'))
        return



