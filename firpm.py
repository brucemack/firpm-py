"""
This file contains a Python translation of the FORTRAN code contained in the 
famous 1973 article entitled "A Computer Program for Designing Linear
Phase Digital Filters" by McClellan, Parks, and Rabiner.  

A PDF copy of the paper is here: https://web.ece.ucsb.edu/Faculty/Rabiner/ece259/Reprints/062_computer%20program.pdf
James McClellan's original Master's thesis from Rice Uiversity is here: https://repository.rice.edu/server/api/core/bitstreams/a924e584-8512-4852-9801-c602985dc0da/content
A helpful reference to the FORTRAN code is here: https://michaelgellis.tripod.com/dsp/pgm21.html

I have tried to keep things as similar to the original FORTRAN as possible
in order to simplify the translation process and to maintain a strong 
cross-reference to the original paper. Therefore, this will not 
necessarily look like "good" Python code.

Bruce MacKinnon (bruce at mackinnon dot com)
"""
import math

PI = 4.0 * math.atan(1.0)
PI2 = 2.0 * PI

# Used to improve consistency with FORTRAN DO loops which 
# iterate across a closed range of [start, end], unlike 
# Python which uses an open range of [start, end)
def CLOSED_RANGE(x, y):
    return range(x, y + 1)

# A 1-based vector to make things easier when porting from FORTRAN 77
class Vector:
    def __init__(self, dim: int = 0):
        self.dim = dim
        self.data = [None] * self.dim
    def fill(self, fill: list):
        self.dim = len(fill)
        self.data = [None] * self.dim
        if not fill is None:            
            for i in range(0, self.dim):
                self.set(i + 1, fill[i])
    def check_index(self, n):
        if n < 1 or n > self.dim:
            raise Exception("Vector index out of bounds")
    def set(self, n, value):
        self.check_index(n)
        self.data[n-1] = value
    def get(self, n):
        self.check_index(n)
        return self.data[n-1]
    def dump(self):
        for i in range(0, self.dim):
            if not self.data[i] is None:
                print(i+1, self.data[i])
    def as_list(self):
        """
        result = []
        for e in self.data:
            # Stop at the first empty element
            if e is None:
                break
            result.append(e)
        return result
        """
        return self.data 
    
class IntVector(Vector):
    def __init__(self, dim):
        Vector.__init__(self, dim)
    def set(self, n, value):
        if not isinstance(value, int):
            raise Exception("Invalid type")
        Vector.set(self, n, value)

# -----------------------------------------------------------------------
# FUNCTION: EFF
#   FUNCTION TO CALCULATE THE DESIRED MAGNITUDE RESPONSE
#   AS A FUNCTION OF FREQUENCY.
#   AN ARBITRARY FUNCTION OF FREQUENCY CAN BE
#   APPROXIMATED IF THE USER REPLACES THIS FUNCTION
#   WITH THE APPROPRIATE CODE TO EVALUATE THE IDEAL
#   MAGNITUDE.  NOTE THAT THE PARAMETER FREQ IS THE
#   VALUE OF NORMALIZED FREQUENCY NEEDED FOR EVALUATION.
# -----------------------------------------------------------------------
def eff(freq, fx: Vector, wtx: Vector, lband, jtype):
    if jtype == 2:
        return fx.get(lband) * freq
    else:
        return fx.get(lband)

# -----------------------------------------------------------------------
# FUNCTION: WATE
#   FUNCTION TO CALCULATE THE WEIGHT FUNCTION AS A FUNCTION
#   OF FREQUENCY.  SIMILAR TO THE FUNCTION EFF, THIS FUNCTION CAN
#   BE REPLACED BY A USER-WRITTEN ROUTINE TO CALCULATE ANY
#   DESIRED WEIGHTING FUNCTION.
# -----------------------------------------------------------------------
def wate(freq, fx: Vector, wtx: Vector, lband, jtype):
    if jtype == 2:
        if fx.get(lband) < 0.0001:
            return wtx.get(lband)
        else:
            return wtx.get(lband) / freq
    else:
        return wtx.get(lband)

# -----------------------------------------------------------------------
# FUNCTION: D
#   FUNCTION TO CALCULATE THE LAGRANGE INTERPOLATION
#   COEFFICIENTS FOR USE IN THE FUNCTION GEE.
# -----------------------------------------------------------------------
def d(k: int, n: int, m: int, x: Vector):
    d = 1.0
    q = x.get(k)
    for l in CLOSED_RANGE(1, m):
        # Step size is m
        for j in range(l, n + 1, m):
            if j - k != 0:
                d = 2.0 * d * (q - x.get(j))
    d = 1.0 / d    
    return d

# -----------------------------------------------------------------------
# FUNCTION: GEE
#   FUNCTION TO EVALUATE THE FREQUENCY RESPONSE USING THE
#   LAGRANGE INTERPOLATION FORMULA IN THE BARYCENTRIC FORM
# -----------------------------------------------------------------------
def gee(k: int, n: int, x: Vector, y: Vector, ad: Vector, grid: Vector):
    """
    Parameters:
        k: Index into the frequency grid to be used for the evaluation.
        n: The number of basis functions being evaluated.
        x: The abscissa vector.  Allows a grid frequency index k to be converted 
           to cos(2 PI F[k])
        y: UNKNOWN
        ad: UNKNOWN
        grid: The grid itself
    """
    p = 0.0
    xf = grid.get(k) 
    xf = math.cos(PI2 * xf)
    d = 0.0 
    for j in CLOSED_RANGE(1, n):
        c = xf - x.get(j)
        c = ad.get(j) / c
        d = d + c 
        p = p + c * y.get(j)
    return p / d 

# -----------------------------------------------------------------------
# SUBROUTINE: REMEZ
#   THIS SUBROUTINE IMPLEMENTS THE REMEZ EXCHANGE ALGORITHM
#   FOR THE WEIGHTED CHEBYSHEV APPROXIMATION OF A CONTINUOUS
#   FUNCTION WITH A SUM OF COSINES.  INPUTS TO THE SUBROUTINE
#   ARE A DENSE GRID WHICH REPLACES THE FREQUENCY AXIS, THE
#   DESIRED FUNCTION ON THIS GRID, THE WEIGHT FUNCTION ON THE
#   GRID, THE NUMBER OF COSINES, AND AN INITIAL GUESS OF THE
#   EXTREMAL FREQUENCIES.  THE PROGRAM MINIMIZES THE CHEBYSHEV
#   ERROR BY DETERMINING THE BEST LOCATION OF THE EXTREMAL
#   FREQUENCIES (POINTS OF MAXIMUM ERROR) AND THEN CALCULATES
#   THE COEFFICIENTS OF THE BEST APPROXIMATION.
#-----------------------------------------------------------------------
def remez(ngrid: int, nfcns: int, grid: Vector, des: Vector, wt: Vector, iext: Vector):
    
    itrmax = int(25)
    devl = -1.0
    nz = int(nfcns + 1)
    nzz = int(nfcns + 2)
    niter = int(0)
    x = Vector(nzz)
    y = Vector(66)
    ad = Vector(nz) 
    comp = 0
    ynz = 0
    a = Vector(65) 
    p = Vector(65) 
    q = Vector(65) 
    k = int(0)
    alpha = Vector(66)

    # Main iteration loop
    # LINE 100 
    while True: 

        iext.set(nzz, ngrid + 1)
        niter = niter + 1

        # Look for iteration max for break out
        if niter > itrmax: 
            break
        
        # This appears to be a performance optimization.  Builds
        # an array of the value of the cosine basis functions
        # sampled at the currently-assumed extremal frequencies.
        for j in CLOSED_RANGE(1, nz):
            jxt = iext.get(j) 
            dtemp = grid.get(jxt) 
            dtemp = math.cos(dtemp * PI2)
            x.set(j, dtemp)
        
        # FORTRAN variables that start with J are implicitly integers
        jet = int((nfcns - 1) / 15 + 1)
        # DO LOOP 120
        for j in CLOSED_RANGE(1, nz):
            ad.set(j, d(j, nz, jet, x))

        dnum = 0.0
        dden = 0.0
        k = int(1)

        for j in CLOSED_RANGE(1, nz):
            l = iext.get(j)
            dtemp = ad.get(j) * des.get(l)
            dnum = dnum + dtemp 
            dtemp = float(k) * ad.get(j) / wt.get(l)
            dden = dden + dtemp
            k = int(-k)        

        dev = dnum / dden 
        nu = 1
        if dev > 0.0: 
            nu = -1 
        dev = -float(nu) * dev
        k = nu 

        # DO LOOP 140
        for j in CLOSED_RANGE(1, nz):
            l = iext.get(j)
            dtemp = float(k) * dev / wt.get(l) 
            y.set(j, des.get(l) + dtemp)
            k = -k 

        if not dev > devl: 
            raise Exception("Error")
        
        # LINE 150
        devl = dev 
        jchnge = 0 
        k1 = iext.get(1)
        knz = iext.get(nz) 
        klow = 0
        nut = -nu
        j = int(1)
        l = int(0)

        # SEARCH FOR THE EXTREMAL FREQUENCIES OF THE BEST
        # APPROXIMATION

        # NOTE: We are using this strange logic structure because
        # of the heavy use of GOTOs in the original FORTRAN code. This
        # is the most straight-forward way I can think of to 
        # replicate the logic without a major re-structuring.

        # Starting point
        GOTO_LINE = 200

        while True:

            if GOTO_LINE == 200:

                if j == nzz:
                    ynz = comp 
                if j >= nzz:
                    GOTO_LINE = 300
                    continue            
                kup = iext.get(j + 1)
                l = iext.get(j) + 1
                nut = -nut 
                if j == 2:
                    y1 = comp 
                comp = dev 

                if l >= kup:
                    GOTO_LINE = 220
                    continue

                err = gee(l, nz, x, y, ad, grid)
                err = (err - des.get(l)) * wt.get(l)
                dtemp = float(nut) * err - comp 
                if dtemp < 0.0: 
                    GOTO_LINE = 220
                    continue
                comp = float(nut) * err 
                GOTO_LINE = 210
                continue 

            elif GOTO_LINE == 210:
                l = l + 1 
                if l >= kup:
                    GOTO_LINE = 215
                    continue
                err = gee(l, nz, x, y, ad, grid)
                err = (err - des.get(l)) * wt.get(l)
                dtemp = float(nut) * err - comp 
                if dtemp < 0.0:
                    GOTO_LINE = 215
                    continue
                comp = float(nut) * err
                GOTO_LINE = 210
                continue

            elif GOTO_LINE == 215:
                iext.set(j, l - 1)
                j = j + 1
                klow = l - 1
                jchnge = jchnge + 1
                GOTO_LINE = 200
                continue

            elif GOTO_LINE == 220:
                l = l - 1
                GOTO_LINE = 225
                continue

            elif GOTO_LINE == 225:
                l = l - 1
                if l <= klow:
                    GOTO_LINE = 250
                    continue 
                err = gee(l, nz, x, y, ad, grid)
                err = (err - des.get(l)) * wt.get(l)
                dtemp = float(nut) * err - comp 
                if dtemp > 0.0:
                    GOTO_LINE = 230
                    continue 
                if jchnge <= 0.0:
                    GOTO_LINE = 225
                    continue
                GOTO_LINE = 260
                continue 

            elif GOTO_LINE == 230:
                comp = float(nut) * err 
                GOTO_LINE = 235
                continue 

            elif GOTO_LINE == 235:
                l = l - 1 
                if l <= klow:
                    GOTO_LINE = 240
                    continue 
                err = gee(l, nz, x, y, ad, grid)
                err = (err - des.get(l)) * wt.get(l)
                dtemp = float(nut) * err - comp 
                if dtemp <= 0.0:
                    GOTO_LINE = 240
                    continue 
                comp = float(nut) * err 
                GOTO_LINE = 235
                continue 

            elif GOTO_LINE == 240:
                klow = iext.get(j)
                iext.set(j, l + 1)
                j = j + 1
                jchnge = jchnge + 1
                GOTO_LINE = 200
                continue 

            elif GOTO_LINE == 250:
                l = iext.get(j) + 1
                if jchnge > 0:
                    GOTO_LINE = 215
                    continue 
                GOTO_LINE = 255
                continue 

            elif GOTO_LINE == 255:
                l = l + 1 
                if l >= kup:
                    GOTO_LINE = 260
                    continue 
                err = gee(l, nz, x, y, ad, grid)
                err = (err - des.get(l)) * wt.get(l)
                dtemp = float(nut) * err - comp 
                if dtemp <= 0.0:
                    GOTO_LINE = 255
                    continue 
                comp = float(nut) * err 
                GOTO_LINE = 210
                continue 

            elif GOTO_LINE == 260:
                klow = iext.get(j) 
                j = j + 1
                GOTO_LINE = 200
                continue 

            elif GOTO_LINE == 300:
                if j > nzz:
                    GOTO_LINE = 320
                    continue 
                if k1 > iext.get(1):
                    k1 = iext.get(1) 
                if knz < iext.get(nz):
                    knz = iext.get(nz) 
                nut1 = nut 
                nut = -nu 
                l = 0
                kup = k1 
                comp = ynz * 1.00001
                luck = 1 
                GOTO_LINE = 310
                continue 

            elif GOTO_LINE == 310:
                l = l + 1
                if l >= kup:
                    GOTO_LINE = 315
                    continue 
                err = gee(l, nz, x, y, ad, grid)
                err = (err - des.get(l)) * wt.get(l)
                dtemp = float(nut) * err - comp 
                if dtemp <= 0.0:
                    GOTO_LINE = 310
                    continue 
                comp = float(nut) * err 
                j = nzz 
                GOTO_LINE = 210
                continue 

            elif GOTO_LINE == 315:
                luck = 6 
                GOTO_LINE = 325
                continue 

            elif GOTO_LINE == 320:
                if luck > 9:
                    GOTO_LINE = 350
                    continue 
                if comp > y1:
                    y1 = comp 
                k1 = iext.get(nzz)
                GOTO_LINE = 325
                continue 

            elif GOTO_LINE == 325:
                l = ngrid + 1 
                klow = knz 
                nut = -nut1 
                comp = y1 * 1.00001 
                GOTO_LINE = 330
                continue 

            elif GOTO_LINE == 330:
                l = l - 1 
                if l <= klow:
                    GOTO_LINE = 340
                    continue 
                err = gee(l, nz, x, y, ad, grid)
                err = (err - des.get(l)) * wt.get(l)
                dtemp = float(nut) * err - comp 
                if dtemp <= 0.0:
                    GOTO_LINE = 330
                    continue 
                j = nzz 
                comp = float(nut) * err 
                luck = luck + 10
                GOTO_LINE = 235
                continue 

            elif GOTO_LINE == 340:
                if luck == 6:
                    GOTO_LINE = 370
                    continue 
                for j in CLOSED_RANGE(1, nfcns):
                    nzzmj = nzz - j 
                    nzmj = nz - j
                    iext.set(nzzmj, iext.get(nzmj))
                iext.set(1, k1)
                GOTO_LINE = 100
                continue 

            elif GOTO_LINE == 350:
                kn = iext.get(nzz)
                # DO LOOP 360
                for j in CLOSED_RANGE(1, nfcns):
                    iext.set(j, iext.get(j + 1))
                iext.set(nz, kn)
                GOTO_LINE = 100
                continue 

            elif GOTO_LINE == 370:
                if jchnge > 0: 
                    GOTO_LINE = 100
                    continue             
                GOTO_LINE = 400
                continue 

            elif GOTO_LINE == 100 or GOTO_LINE == 400:
                # Here is where we can break out of this strange contstruct
                break

            else:
                raise Exception("Invalid GOTO_LINE")

        # Decide how to handle the master iteration based on the outcome of the
        # previous GOTO construct            
        if GOTO_LINE == 100:
            continue 
        elif GOTO_LINE == 400:
            break
        else:
            raise Exception("Invalid GOTO_LINE")

    # CALCULATION OF THE COEFFICIENTS OF THE BEST APPROXIMATION
    # USING THE INVERSE DISCRETE FOURIER TRANSFORM
    
    # LINE 400

    nm1 = nfcns - 1 
    fsh = 1.0e-6
    gtemp = grid.get(1) 
    x.set(nzz, -2.0) 
    cn = 2 * nfcns - 1
    delf = 1.0 / cn 
    l = 1 
    kkk = 0
    if grid.get(1) < 0.01 and grid.get(ngrid) > 0.49:
        kkk = 1
    if nfcns <= 3:
        kkk = 1 
    if not kkk == 1:
        dtemp = math.cos(PI2 * grid.get(1))
        dnum = math.cos(PI2 * grid.get(ngrid))
        aa = 2.0 / (dtemp - dnum)
        bb = -(dtemp + dnum) / (dtemp - dnum)

    # Line 405 
    # DO LOOP 430
    for j in CLOSED_RANGE(1, nfcns):
        ft = j - 1 
        ft = ft * delf
        xt = math.cos(PI2 * ft)
        if not kkk == 1: 
            xt = (xt - bb) / aa 
            xt1 = math.sqrt(1.0 - xt * xt)
            ft = math.atan2(xt1, xt) / PI2

        # There is another set of complicated GOTOs in the original code
        # adopting the strange GOTO_LINE construct to replicate.
        GOTO_LINE = 410

        while True:

            if GOTO_LINE == 410:
                xe = x.get(l) 
                if xt > xe:
                    GOTO_LINE = 420
                    continue 
                if (xe - xt) < fsh:
                    GOTO_LINE = 415 
                    continue 
                l = l + 1 
                GOTO_LINE = 410 
                continue

            elif GOTO_LINE == 415:
                a.set(j, y.get(l))
                GOTO_LINE = 425 
                continue 

            elif GOTO_LINE == 420:
                if (xt - xe) < fsh:
                    GOTO_LINE = 415 
                    continue 
                grid.set(1, ft)
                a.set(j, gee(1, nz, x, y, ad, grid))
                GOTO_LINE = 425 
                continue 

            elif GOTO_LINE == 425:
                if l > 1:
                    l = l - 1
                GOTO_LINE = 430 
                continue

            elif GOTO_LINE == 430:
                # The only way to get out 
                break 
            
            else:
                raise Exception("Invalid GOTO_LINE")
    # Line 430 
    grid.set(1, gtemp)
    dden = PI2 / cn 

    # DO LOOP 510
    for j in CLOSED_RANGE(1, nfcns):
        dtemp = 0.0 
        dnum = j - 1
        dnum = dnum * dden 
        if not nm1 < 1: 
            # DO LOOP 500
            for k in CLOSED_RANGE(1, nm1):
                dak = a.get(k + 1)
                dk = k 
                dtemp = dtemp + dak * math.cos(dnum * dk)
        dtemp = 2.0 * dtemp + a.get(1)
        alpha.set(j, dtemp)

    # DO LOOP 550          
    for j in CLOSED_RANGE(2, nfcns):
        alpha.set(j, 2.0 * alpha.get(j) / cn)

    # Close to line 550
    alpha.set(1, alpha.get(1) / cn )

    # Or else jump to line 545
    if not kkk == 1: 
        p.set(1, 2.0 * alpha.get(nfcns) * bb + alpha.get(nm1))
        p.set(2, 2.0 * aa * alpha.get(nfcns))
        q.set(1, alpha.get(nfcns - 2) - alpha.get(nfcns))

        # DO LOOP 540 
        for j in CLOSED_RANGE(2, nm1):

            if not j < nm1:
                aa = 0.5 * aa 
                bb = 0.5 * bb 
            p.set(j + 1, 0.0)

            # DO LOOP 520 
            for k in CLOSED_RANGE(1, j):
                a.set(k, p.get(k))
                p.set(k, 2.0 * bb * a.get(k))
            
            p.set(2, p.get(2) + a.get(1) * 2.0 * aa)
            jm1 = j - 1 

            # DO LOOP 525 
            for k in CLOSED_RANGE(1, jm1):
                p.set(k, p.get(k) + q.get(k) + aa * a.get(k + 1))
            
            jp1 = j + 1 

            # DO LOOP 530
            for k in CLOSED_RANGE(3, jp1):
                p.set(k, p.get(k) + aa * a.get(k - 1))

            if not j == nm1:

                # DO LOOP 535
                for k in CLOSED_RANGE(1, j): 
                    q.set(k, -a.get(k))
                
                nf1j = nfcns - 1 - j

                q.set(1, q.get(1) + alpha.get(nf1j))

        # Line 540
        # DO LOOP 543
        for j in CLOSED_RANGE(1, nfcns):
            alpha.set(j, p.get(j))
            
    # LINE 545
    if not nfcns > 3:
        alpha.set(nfcns + 1, 0.0)
        alpha.set(nfcns + 2, 0.0)

    # TODO: BUILD THE FINAL G(F) FOR RETURN
    return dev, alpha

def design(nfilt: int, jtype: int, nbands: int, edges: list, gains: list, weights: list, lgrid: int = 16):
    """

    General Notes:
    
    Frequencies are specified in normalized form [0:1] where the frequency 
    is interpreted as cycles per sample.  The Nyquist frequency is at 0.5
    cycles/sample.

    Parameters:
        nfilt: Filter length (taps)
        jtype: 
        nbands: Numer of bands
        edges: Pairs of frequencies defining the start and end point of each 
            band. Maximum of 10 bands allowed.
        gains: Response per band.  One entry per band.
        weights: Weight per band.  One entry per band.
        lgrid: The grid density used for estimation.

    Returns:
        The impulse response (list of nfilt coefficients)
        The deviation achieved
    """

    edge = Vector()
    edge.fill(edges)
    fx = Vector()
    fx.fill(gains)
    wtx = Vector()
    wtx.fill(weights)

    nfmax:int = int(128)

    # Validation
    if not (nfilt <= nfmax and nfilt >= 3):
        raise Exception("Invalid nfilt specified")

    # LINE 115
    jb = 2 * nbands
    # LINE 120
    if not (jtype > 0 and jtype <= 3):
        raise Exception("Invalid jtype specified")
    # LINE 125
    neg = int(1)
    # TODO: DOC
    # Multiple pass-band/stop-band filters need to antisymmetric
    if jtype == 1: 
        neg = 0
    # This will be zero if the filter length is even
    # and one if its odd.
    nodd = int(nfilt / 2)
    nodd = int(nfilt - 2 * nodd)
    # nfcns is the number of basis functions used to approximate the 
    # magnitude response.
    nfcns = int(nfilt / 2)
    # TODO: Document this
    if nodd == 1 and neg == 0: 
        nfcns = nfcns + 1

    # Setup the dense grid of frequencies.
    grid = Vector((nfilt + 1) * int(lgrid / 2) + 1)
    des = Vector((nfilt + 1) * int(lgrid / 2))
    wt = Vector((nfilt + 1) * int(lgrid / 2))
    # Start off the grid at the lower-boundary of the first band
    grid.set(1, edge.get(1))
    # delf is the delta in frequency between grid points.
    delf = lgrid * nfcns
    delf = 0.5 / delf
    if not neg == 0:
        if edge.get(1) < delf:
            grid.set(1, delf)
    
    # LINE 135
    j = int(1)
    l = int(1)
    lband = int(1)
    iext = IntVector(1045)
    # TODO: CHECK THIS
    h = Vector(round(nfilt / 2) + 1)

    # This is the iteration across the bands.  The index "l"
    # points to the current band.
    while True:

        # LINE 140
        # Upper bound of band we are working on
        fup = edge.get(l + 1)

        # This is the iteration across the grid within a single band.
        # The "j" index points to the current grid point.
        while True:
            # Line 145
            temp = grid.get(j) 

            # CALCULATE THE DESIRED MAGNITUDE RESPONSE AND THE WEIGHT
            # FUNCTION ON THE GRID
            des.set(j, eff(temp, fx, wtx, lband, jtype))
            wt.set(j, wate(temp, fx, wtx, lband, jtype))

            j = j + 1
            grid.set(j, temp + delf)

            # Check to see if we've over-stepped the bands boundary.
            # If so break out.
            if grid.get(j) > fup:
                break
        
        # LINE 150
        # Back-up and correct the grid point based on the upper frequency
        # of the band.
        grid.set(j - 1,  fup)
        des.set(j - 1, eff(fup, fx, wtx, lband, jtype))
        wt.set(j - 1, wate(fup, fx, wtx, lband, jtype))
        # Advance to the next band
        lband = lband + 1
        l = l + 2

        # Look for exit from bands loop
        if lband > nbands: 
            break

        # The first point in the band's grid is the lower-boundary 
        # of the band.
        grid.set(j, edge.get(l))

    # Line 160
    ngrid = j - 1 
    if not neg != nodd:
        # TODO: Document what is happening here.
        if grid.get(ngrid) > (0.5 - delf):
            ngrid = ngrid - 1

    # SET UP A NEW APPROXIMATION PROBLEM WHICH IS EQUIVALENT
    # TO THE ORIGINAL PROBLEM.
    #
    # Adjust the problem based on the type of filter and
    # symmetry.
    if neg == 0:
        if not nodd == 1:
            for j in CLOSED_RANGE(1, ngrid):
                change = math.cos(PI * grid.get(j))
                des.set(j, des.get(j) / change)
                wt.set(j, wt.get(j) * change)
    else:
        if not nodd == 1:
            for j in CLOSED_RANGE(1, ngrid):
                change = math.sin(PI * grid.get(j))
                des.set(j, des.get(j) / change)
                wt.set(j, wt.get(j) * change)
        else:
            for j in CLOSED_RANGE(1, ngrid):
                change = math.sin(PI2 * grid.get(j))
                des.set(j, des.get(j) / change)
                wt.set(j, wt.get(j) * change)

    # INITIAL GUESS FOR THE EXTREMAL FREQUENCIES--EQUALLY
    # SPACED ALONG THE GRID
    # LINE 200
    temp = float(ngrid - 1) / float(nfcns)
    for j in CLOSED_RANGE(1, nfcns):
        xt = j - 1
        # TODO: Need to review whether the conversion to integer
        # is correct here.
        iext.set(j, int(xt * temp + 1.0))
    iext.set(nfcns + 1, ngrid)
    nm1 = nfcns - 1
    nz = nfcns + 1

    # Call the big function
    dev, alpha = remez(ngrid, nfcns, grid, des, wt, iext)

    # Implement equations (9) - (12) 

    if neg == 0:
        # Line 300
        if not nodd == 0: 
            # DO LOOP 305
            for j in CLOSED_RANGE(1, nm1):
                nzmj = nz - j
                h.set(j, 0.5 * alpha.get(nzmj))
            h.set(nfcns, alpha.get(1))
        else:
            # Line 310
            h.set(1, 0.25 * alpha.get(nfcns))
            # DO LOOP 315
            for j in CLOSED_RANGE(2, nm1):
                nzmj = nz - j
                nf2j = nfcns + 2 - j
                h.set(j, 0.25 * (alpha.get(nzmj) + alpha.get(nf2j)))
            h.set(nfcns, 0.5 * alpha.get(1) + 0.25 * alpha.get(2))
    # Line 320
    else:
        if not nodd == 0:
            h.set(1, 0.25 * alpha.get(nfcns))
            h.set(2, 0.25 * alpha.get(nm1))
            # DO LOOP 325
            for j in CLOSED_RANGE(3, nm1):
                nzmj = nz - j
                nf3j = nfcns + 3 - j
                h.set(j, 0.25 * (alpha.get(nzmj) - alpha.get(nf3j)))
            h.set(nfcns, 0.5 * alpha.get(1) - 0.25 * alpha.get(3))
            h.set(nz, 0)            
        else:
            # Line 330
            h.set(1, 0.25 * alpha.get(nfcns))
            # DO LOOP 335
            for j in CLOSED_RANGE(2, nm1):
                nzmj = nz - j
                nf2j = nfcns + 2 - j
                h.set(j, 0.25 * (alpha.get(nzmj) - alpha.get(nf2j)))
            h.set(nfcns, 0.5 * alpha.get(1) - 0.25 * alpha.get(2))         

    # Make a complete impulse response 
    impulse = []
    half_size = int(nfilt / 2)
    # Left half
    for i in CLOSED_RANGE(1, half_size):
        impulse.append(h.get(i))
    # Center
    if nodd == 1:
        impulse.append(h.get(half_size + 1))
    # Right half
    for i in CLOSED_RANGE(1, half_size):
        tap = h.get(half_size - i + 1)
        # TODO: CHECK THIS LOGIC
        if jtype == 2 or jtype == 3:
            tap = -1 * tap 
        impulse.append(tap)

    return impulse, dev
