from INGRID import ingrid
#from Interpol.Setup_Grid_Data import *
import e2dgrid as e
from  matplotlib.pyplot import figure,ion
#import numpy as np
#from copy import copy
ion()


def set_gridue_manual(ingridobj,rm,zm,nx,ny,ixpt1,ixpt2,iysptrx):
        from numpy import zeros,sqrt
        """
        Prepares a ``gridue_settings`` dictionary with required data
        for writing a gridue file.

        Parameters
        ----------

        Returns
        -------

        """

        # RECALL: ingridobj.rm has FORTRAN style ordering (columns are accessed via the first entry)
        # Getting relevant values for gridue file
        ixrb = len(rm) - 2
        nxm = len(rm) - 2
        nym = len(rm[0]) - 2


        psi = zeros((nxm + 2, nym + 2, 5), order='F')
        br = zeros((nxm + 2, nym + 2, 5), order='F')
        bz = zeros((nxm + 2, nym + 2, 5), order='F')
        bpol = zeros((nxm + 2, nym + 2, 5), order='F')
        bphi = zeros((nxm + 2, nym + 2, 5), order='F')
        b = zeros((nxm + 2, nym + 2, 5), order='F')

        rb_prod = ingridobj.PsiUNorm.rcenter * ingridobj.PsiUNorm.bcenter

        for i in range(len(b)):
            for j in range(len(b[0])):
                for k in range(5):
                    _r = rm[i][j][k]
                    _z = zm[i][j][k]

                    _psi = ingridobj.PsiUNorm.get_psi(_r, _z)
                    _br = ingridobj.PsiUNorm.get_psi(_r, _z, tag='vz') / _r
                    _bz = -ingridobj.PsiUNorm.get_psi(_r, _z, tag='vr') / _r
                    _bpol = sqrt(_br ** 2 + _bz ** 2)
                    _bphi = rb_prod / _r
                    _b = sqrt(_bpol ** 2 + _bphi ** 2)

                    psi[i][j][k] = _psi
                    br[i][j][k] = _br
                    bz[i][j][k] = _bz
                    bpol[i][j][k] = _bpol
                    bphi[i][j][k] = _bphi
                    b[i][j][k] = _b

        ingridobj.CurrentTopology.gridue_settings = {
            'nxm': nxm, 'nym': nym, 'ixpt1': ixpt1, 'ixpt2': ixpt2, 'iyseptrx1': iysptrx,
            'rm': rm, 'zm': zm, 'psi': psi, 'br': br, 'bz': bz, 'bpol': bpol, 'bphi': bphi, 'b': b
        }


i=ingrid.Ingrid()
i.settings['grid_settings']['num_xpt'] = 1
i.settings['eqdsk'] = '/Users/holma2/Dropbox (Aalto)/UEDGE/INGRID/INGRID/INGRID/seq#1/eqdsk'
i.settings['grid_settings']['rxpt'] = 0
i.settings['grid_settings']['zxpt'] = -1
i.StartSetup()
i.SetTopology('LSN')
set_gridue_manual(i,*e.grid('/Users/holma2/Dropbox (Aalto)/UEDGE/E2D-EIR_integration/exporte2dgrid/aholm_mar1719_11'))
i.ExportGridue()

# Create a test case to display the UEDGE solution
from uedge import *
import uedge.contrib.holm10.plot as p
from matplotlib.pyplot import subplots


bbb.gengrid=0
bbb.ngrid=1
bbb.mhdgeo=1
com.geometry='snull'
bbb.ftol=1e20
#bbb.issfon=0
exmain()

# Display INGRID gridue
titles = ['Tot. B-field', 'Psi', 'B_r', 'B_z', 'B_pol', 'B_phi']
params = [com.b, com.psi, com.br, com.bz, com.bpol, com.bphi]
subtitle = ['cell-center', 'SE', 'SW', 'NE', 'NW']

for v in range(5):
    f, ax = subplots(1, 5, figsize=(20,5))
    f.suptitle(titles[v])
    for i in range(5):
        p.heatmap(params[v],s=i,title=subtitle[i],ax=ax[i], zoom='device')


