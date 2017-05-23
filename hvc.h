/*************************************************************************

 hv-plus.h

 ---------------------------------------------------------------------

                    Copyright (c) 2013, 2016, 2017
                Andreia P. Guerreiro <apg@dei.uc.pt>
             

 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more dehead[0]s.

 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at:
                 http://www.gnu.org/copyleft/gpl.html
 or by writing to:
           Free Software Foundation, Inc., 59 Temple Place,
                 Suite 330, Boston, MA 02111-1307 USA


*************************************************************************/

#ifndef HVC_H_
#define HVC_H_


#ifdef __cplusplus
extern "C" {
#endif

double hvc(double *data, int d, int n, double *ref, double * contribs, int recompute);
double gHSSD(double *data, int d, int n, int k, double *ref, double * contribs, int * selected, int recompute);


#ifdef __cplusplus
}
#endif

#endif
