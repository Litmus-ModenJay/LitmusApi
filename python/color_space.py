import math

# Color Vector Conversion & Classification
class CVC():

# RGB, rgb, hexa
    @staticmethod
    def hexa_RGB(hexa):
        return tuple(int(hexa[2*i+1:2*i+3], 16) for i in range(0,3))

    @staticmethod
    def hexa_rgb(hexa):
        return tuple(float(int(hexa[2*i+1:2*i+3], 16))/255 for i in range(0,3))
    
    @staticmethod
    def rgb_RGB(rgb):
        RGB = [-1, -1, -1]
        for i in range(0,3):
            if rgb[i] >= 0 and rgb[i] <= 1 :
                RGB[i] = int(round(rgb[i]*255))
        return tuple(RGB[i] for i in range(0,3))
    
    @staticmethod
    def RGB_rgb(RGB):
        return tuple(float(RGB[i])/255 for i in range(0,3))
    
    @staticmethod
    def RGB_hexa(RGB):
        if RGB[0] < 0 or RGB[1] < 0 or RGB[2] < 0:
            pre =  '@'
        else:
            pre = '#'
        hexa =  pre + "{0}{1}{2}".format("%02X" % (RGB[0]),"%02X" % (RGB[1]),"%02X" % (RGB[2]))
        return hexa.upper()

    @staticmethod
    def rgb_hexa(rgb):
        RGB = rgb_RGB(rgb)
        return RGB_hexa(RGB)

    @staticmethod
    def rgb_parameters(rgb):
        summation = sum(rgb)
        maximum = max(rgb)
        minimum = min(rgb)
        sigma = maximum + minimum
        delta = maximum - minimum
        return (summation, maximum, minimum, sigma, delta)

# Cell
    @staticmethod
    def RGB_Cell(RGB):
        Ch = tuple(int(RGB[i]/32) + 1 for i in range(0,3))
        cell = '{0}{1}{2}'.format(Ch[0],Ch[1],Ch[2])
        return cell

# HSLrgb, HSLV, CMYK
    @staticmethod
    def rgb_HSLrgb(rgb):
        r, g, b = rgb[0], rgb[1], rgb[2]
        x = (2*r-g-b)/2.0
        y = math.sqrt(3.0)*(g-b)/2.0
        z = (r + g + b)/3.0
        if abs(x)<0.00001 and abs(y)<0.00001 :
            hue = 0
        elif y<0 :
            hue = math.degrees(math.atan2(y,x)) + 360
        else :
            hue = math.degrees(math.atan2(y,x))
        chroma = math.sqrt(x*x + y*y)
        if z<0.00001 or z>0.99999 :
            saturation = 0
        else : 
            saturation = chroma / z
        return (hue, saturation, z, chroma)
    
    @staticmethod
    def rgb_HSLV(rgb):
        r, g, b = rgb[0], rgb[1], rgb[2]
        maximum = max(rgb)
        minimum = min(rgb)
        sigma = maximum + minimum
        delta = maximum - minimum
        light = sigma/2.0
        if delta == 0 :
            hue, sl, sv = 0, 0, 0
        else :
            if maximum ==  r :
                hue = (float((g-b)/delta) % 6) * 60
            elif maximum == g :
                hue = (float((b-r)/delta) % 6 + 2) * 60
            else :
                hue = (float((r-g)/delta) % 6 + 4) * 60
            sl = delta / (1 - abs(2*light - 1))
            if maximum :
                sv = delta/maximum
            else:
                sv = 0
        return (hue, sl, light, sv, maximum, delta)

    @staticmethod
    def rgb_CMYK(rgb):
        r, g, b = rgb[0], rgb[1], rgb[2]
        k = 1 - max(rgb)
        if k == 1 :
            c, m, y = 0, 0, 0
        else :
            c = (1 - r - k) / (1 - k)
            m = (1 - g - k) / (1 - k)
            y = (1 - b - k) / (1 - k)
        return (c, m, y, k)

# XYZ
    @staticmethod
    def rgb_XYZ(rgb, profile):
        rgb_linearized = [0,0,0]
        for i in range(0,3):
            if rgb[i] < 0.04045:
                rgb_linearized[i] = rgb[i] / 12.92
            else:
                rgb_linearized[i] = ((rgb[i] + 0.055) / 1.055) ** 2.4
        r, g, b = rgb_linearized[0], rgb_linearized[1], rgb_linearized[2]
        if profile == "sRGB" :
            X = 0.4124564*r + 0.3575761*g + 0.1804375*b
            Y = 0.2126729*r + 0.7151522*g + 0.0721750*b
            Z = 0.0193339*r + 0.1191920*g + 0.9503041*b
        elif profile == "AdobeRGB" :
            X = 0.5767309*r + 0.1855540*g + 0.1881852*b
            Y = 0.2973769*r + 0.6273491*g + 0.0752741*b
            Z = 0.0270343*r + 0.0706872*g + 0.9911085*b
        elif profile == "CIERGB" :
            X = 0.4887180*r + 0.3106803*g + 0.2006017*b
            Y = 0.1762044*r + 0.8129847*g + 0.0108109*b
            Z = 0.0000000*r + 0.0102048*g + 0.9897952*b
        sum = X + Y + Z
        if sum == 0 :
            x, y = 0.0, 0.0
        else :
            x, y = X / sum, Y / sum
        return (X, Y, Z, sum/3, x, y)                 

    @staticmethod
    def XYZ_rgb(XYZ):
        X, Y, Z = XYZ[0], XYZ[1], XYZ[2]
        r = 3.2404542*X + (-1.5371385)*Y + (-0.4985314)*Z
        g = (-0.9692660)*X + 1.8760108*Y + 0.0415560*Z
        b = 0.0556434*X + (-0.2040259)*Y + 1.0572252*Z
        rgb = (r,g,b)
        rgb_dl = [0,0,0]
        for i in range(0,3):
            if rgb[i] <= 0.0031308:
                rgb_dl[i] = rgb[i] * 12.92     
            else:
                rgb_dl[i] = (rgb[i])**(1/2.4)*1.055 - 0.055
            if abs(rgb_dl[i]) < 0.000001:
                rgb_dl[i] = 0.0
            if abs(rgb_dl[i]-1.0) < 0.000001:
                rgb_dl[i] = 1.0 
    
        return (rgb_dl[0],rgb_dl[1],rgb_dl[2])


# F fucntion for Lab Luv
    @staticmethod
    def F(value):
        delta = 6/29
        if value > delta**3.0:
            transform = value**(1.0/3.0)
        else:
            transform = value / 3 / (delta**2.0) + 4 / 29
        return transform
    
    @staticmethod
    def F_inv(value):
        delta = 6/29
        if value > delta:
            transform = value**3.0
        else:
            transform = (value - 4/29) * 3 * (delta**2.0)
        return transform

    @staticmethod
    def XYZ_Labuv(XYZ, illuminant):
        X, Y, Z = XYZ[0], XYZ[1], XYZ[2]
        if illuminant == "D50_2" :
            Xn, Yn, Zn = 0.966797, 1.0000, 0.825188
        elif illuminant == "D65_2" :
            Xn, Yn, Zn = 0.95047, 1.0000, 1.08883
        elif illuminant == "E" :
            Xn, Yn, Zn = 1.0000, 1.0000, 1.0000
        
        xn, yn = Xn / (Xn + Yn + Zn), Yn / (Xn + Yn + Zn)
        un, vn = 4*xn / ((-2)*xn + 12*yn + 3), 9*yn / ((-2)*xn + 12*yn + 3)
        if X+Y+Z == 0:
            L = 0.0001
            a, b = 0, 0
            u, v = 0, 0
            ua, va = 0, 0
        else :
            L = (116 * CVC.F(Y/Yn) - 16) / 100
            a = (500 * (CVC.F(X/Xn) - CVC.F(Y/Yn))) / 100
            b = (200 * (CVC.F(Y/Yn) - CVC.F(Z/Zn))) / 100
            u = 4 * X / (X + 15*Y + 3*Z)
            v = 9 * Y / (X + 15*Y + 3*Z)
            ua, va = 13*L*(u-un), 13*L*(v-vn)
        if abs(a)<0.00001 and abs(b)<0.00001 :
            H_lab = 0
        elif b<0 :
            H_lab = math.degrees(math.atan2(b,a)) + 360
        else :
            H_lab = math.degrees(math.atan2(b,a))
        if abs(ua)<0.00001 and abs(va)<0.00001 :
            H_luv = 0
        elif va<0 :
            H_luv = math.degrees(math.atan2(va,ua)) + 360
        else :
            H_luv = math.degrees(math.atan2(va,ua))
        C_lab = math.sqrt(a*a + b*b)
        C_luv = math.sqrt(ua*ua + va*va)
        if L<0.00001 or L>0.99999 :
            S_lab = 0
            S_luv = 0
            if L>0.99999:
                L = 0.99999
        else :
            S_lab = C_lab / L
            S_luv = C_luv / L
        return (L, a, b, H_lab, S_lab, C_lab, ua, va, H_luv, S_luv, C_luv)
    
# Lab coordinate
    @staticmethod
    def XYZ_Lab(XYZ):
        X, Y, Z = XYZ[0], XYZ[1], XYZ[2]
        Xn, Yn, Zn = 0.95043, 1.0000, 1.08889

        if X+Y+Z == 0:
            L = 0.0000001
            a, b = 0.0, 0.0
        else :
            L = (116 * F(Y/Yn) - 16)
            a = (500 * (F(X/Xn) - F(Y/Yn)))
            b = (200 * (F(Y/Yn) - F(Z/Zn)))
        if abs(a)<0.001 and abs(b)<0.001 :
            H_lab = 0.0
        elif b<0 :
            H_lab = math.degrees(math.atan2(b,a)) + 360
        else :
            H_lab = math.degrees(math.atan2(b,a))
        C_lab = math.sqrt(a*a + b*b)
        if L < 0.00001 or L > 99.99999 :
            S_lab = 0.0
            if L > 99.99999:
                L = 99.99999
        else :
            S_lab = C_lab / L
        return (L, a, b, H_lab, C_lab)

    @staticmethod
    def LCH_Lab(LCH):
        a = LCH[1]*math.cos(math.radians(LCH[2]))
        b = LCH[1]*math.sin(math.radians(LCH[2]))
        return (LCH[0], a, b, LCH[1], LCH[2])

    @staticmethod
    def Lab_XYZ(Lab):
        L, a, b = Lab[0], Lab[1], Lab[2]
        Xn, Yn, Zn = 0.95043, 1.0000, 1.08889

        if L == 0:
            X, Y, Z = 0.0, 0.0, 0.0
        else :
            Y = Yn * F_inv((L + 16)/116)
            X = Xn * F_inv((L + 16)/116 + a/500)
            Z = Zn * F_inv((L + 16)/116 - b/200)
        return (X, Y, Z)


# Ypq
    @staticmethod
    def rgb_Ypq(rgb):
        Y = 0.299*rgb[0] + 0.587*rgb[1] + 0.114*rgb[2]
        p = rgb[0]-0.5*rgb[1]-0.5*rgb[2]
        q = 0.8660254*rgb[1]-0.8660254*rgb[2]
        if abs(p) < 0.001 and abs(q) < 0.001 :
            H = 0.0
        elif q < 0 :
            H = math.degrees(math.atan2(q,p)) + 360
        else :
            H = math.degrees(math.atan2(q,p))
        C = math.sqrt(p*p + q*q)
        return (Y, p, q, H, C)

    @staticmethod
    def Ypq_rgb(Ypq):
        r = Ypq[0] + 0.701*Ypq[1] - 0.2730867*Ypq[2]
        g = Ypq[0] - 0.299*Ypq[1] + 0.3042636*Ypq[2]
        b = Ypq[0] - 0.299*Ypq[1] - 0.850437*Ypq[2]
        return (r,g,b)
    
    @staticmethod
    def YCH_Ypg(YCH):
        p = YCH[1]*math.cos(math.radians(YCH[2]))
        q = YCH[1]*math.sin(math.radians(YCH[2]))
        return (YCH[0], p, q, YCH[2], YCH[1])

    @staticmethod
    def YCH_hexa(YCH):
        return rgb_hexa(Ypq_rgb(YCH_Ypg(YCH)))

    @staticmethod
    def hexa_Ypg(hexa):
        return rgb_Ypq(hexa_rgb(hexa))

    @staticmethod
    def Ypq_hexa(Ypq):
        return rgb_hexa(Ypq_rgb(Ypq))

    @staticmethod
    def Ypq_Section(Ypq):
        Href = [9.06,41.89,61.98,95.81,134.92,171.83,188.83,213.75,248.16,283.61,313.71,338.87,369.06]
        Lref = [0.0000,0.2722,0.50196,0.6464,0.75294,0.9287,1.000]
        Cref = [0.0,0.0306,0.0622,0.1706,0.3630,0.5969,0.8313,1.0001]
        
        H, L, C = Ypq[3],Ypq[0],Ypq[4]
        for l in range(len(Lref)-1):
            L0, L1 = Lref[l], Lref[l+1]
            if C < Cref[1] and L >= L0 and L < L1:
                section = 'HacC0Y{0}'.format(l)
            for c in range(1,len(Cref)-1):
                for h in range(len(Href)-1):
                    if H < Href[0]:
                        H = H + 360
                    H0, H1 = Href[h], Href[h+1]
                    L0, L1 = Lref[l], Lref[l+1]
                    C0, C1 = Cref[c], Cref[c+1]         
                    if H >= H0 and H < H1 and L >= L0 and L < L1 and C >= C0 and C < C1:
                        h_new = h+1
                        if h_new == len(Href)-1:
                            h_new = 0
                        section = 'H{0:02d}C{1}Y{2}'.format(h_new,c,l)
        return section

    @staticmethod
    def Section_Wheel(section):
        WHref = ['Red','Orange','Yellow','Lime','Green','Sea','Cyan','Sky',
                'Blue','Violet','Purple','Rose']
        Dref = [['Black','Dark Gray','Gray','Gray','Light Gray','White'],
                ['Black Dark','Dark Gray','Gray','Gray','Light Gray','White Pale'],
                ['Dark','Dark','Dull','Soft','Pale','Pale'],
                ['Dark','Dark','Dull','Soft','Pale','Pale'],
                ['Deep','Deep','Strong','Strong','Light','Light'],
                ['Deep','Deep','Strong','Strong','Light','Light'],
                ['Vivid','Vivid','Vivid','Vivid','Vivid','Vivid']]

        h = section[1:3]
        c = int(section[4])
        l = int(section[6])
        if h == 'ac':
            wheel = 'Achromatic'
            depth = Dref[0][l]
        else:
            wheel = WHref[int(h)]
            depth = Dref[c][l]
        mode = 'Dark'
        if l >= 3:
            mode = 'Light'
        if l == 2 and c == 6 and (int(h) >=  4 and int(h) <= 6):
            mode = 'Dark'
        
        return {'Depth':depth,'Wheel':wheel, 'Mode':mode}

    @staticmethod
    def Section_Group(section):
        Supernova = ['Empty','Red','Orange','Yellow','Green','Blue','Purple',
                    'Pink','Brown','Gray','White','Black']
        Agroup = [11,9,9,9,9,10]
        Cgroup = [[[11,9,9,9,7,10],[8,8,8,7,7,7],[8,8,7,7,7,0],[1,1,7,7,7,0],[1,1,1,0,0,0],[1,1,0,0,0,0]],
                [[11,9,9,9,9,10],[8,8,8,8,2,2],[8,8,8,8,2,2],[8,8,2,2,2,0],[8,2,2,2,2,0],[0,2,2,2,0,0]],
                [[11,9,9,9,9,10],[8,8,8,8,3,3],[8,8,8,8,3,3],[0,8,8,3,3,3],[0,8,3,3,3,3],[0,0,0,3,3,0]],
                [[11,9,9,9,9,10],[4,4,4,4,4,4],[4,4,4,4,4,4],[0,4,4,4,4,4],[0,4,4,4,4,0],[0,0,0,4,4,0]],
                [[11,9,9,9,9,10],[4,4,4,4,4,4],[4,4,4,4,4,4],[4,4,4,4,4,0],[0,4,4,4,4,0],[0,4,4,4,0,0]],
                [[11,9,9,9,9,10],[4,4,4,4,4,4],[4,4,4,4,4,4],[4,4,4,4,4,0],[0,4,4,4,4,0],[0,0,4,4,0,0]],
                [[11,9,9,9,9,10],[4,4,5,5,5,5],[4,4,5,5,5,5],[4,4,5,5,5,0],[0,4,5,5,5,0],[0,0,5,5,0,0]],
                [[11,9,9,9,9,10],[5,5,5,5,5,5],[5,5,5,5,5,0],[5,5,5,5,5,0],[5,5,5,5,0,0],[0,5,5,0,0,0]],
                [[11,9,9,9,9,10],[5,5,5,5,5,5],[5,5,5,5,5,0],[5,5,5,5,0,0],[5,5,5,0,0,0],[5,5,0,0,0,0]],
                [[11,9,9,9,9,10],[6,6,6,6,6,6],[6,6,6,6,6,0],[6,6,6,6,0,0],[6,6,6,0,0,0],[6,6,0,0,0,0]],
                [[11,9,9,9,9,10],[6,6,6,6,6,6],[6,6,6,6,6,0],[6,6,6,6,6,0],[6,6,6,6,0,0],[0,6,6,0,0,0]],
                [[11,9,9,9,7,10],[6,6,6,7,7,7],[6,6,6,7,7,0],[1,6,7,7,0,0],[1,1,7,0,0,0],[0,7,0,0,0,0]]]
        h = (section[1:3])
        l = int(section[6])
        c = int(section[4])
        if h == 'ac':
            group = Agroup[l]
        else:
            group = Cgroup[int(h)][c-1][l]
        
        return Supernova[group]

# GEO coordinate  
    @staticmethod
    def rgb_GEOrgb(rgb):
        HSLrgb= CVC.rgb_HSLrgb(rgb)
        hue = HSLrgb[0]
        z = HSLrgb[2] * 2 - 1
        c = HSLrgb[3]
        radius = (z**2 + c**2)**0.5
        theta = z / radius
        if hue >= 180 :
            longitude = hue - 360
        else:
            longitude = hue
        latitude = math.degrees(math.asin(z))
        return (latitude, longitude, radius)
    
    @staticmethod
    def rgb_GEOHSL(rgb):
        HSLV= CVC.rgb_HSLV(rgb)
        hue = HSLV[0]
        z = HSLV[2] * 2 -1
        c = HSLV[5]
        radius = (z**2 + c**2)**0.5
        theta = z / radius
        if hue >= 180 :
            longitude = hue - 360
        else:
            longitude = hue
        latitude = math.degrees(math.asin(z)) * 0.99
        return (latitude, longitude, radius)

    @staticmethod
    def Labuv_GeoLab(Labuv):
        hue = Labuv[3]
        z = Labuv[0]*2 -1
        c = Labuv[5]
        radius = (z**2 + c**2)**0.5
        theta = z / radius
        if hue >= 180 :
            longitude = hue - 360
        else:
            longitude = hue
        latitude = math.degrees(math.asin(z))

        # radius = ((Labuv[0]-1/2)**2 + Labuv[1]**2 + Labuv[2]**2)**(1.0/2.0)
        return (latitude, longitude, radius)
    
    @staticmethod
    def Labuv_GeoLuv(Labuv):
        hue = Labuv[8]
        z = Labuv[0]*2 -1
        c = Labuv[10]
        radius = (z**2 + c**2)**0.5
        theta = z / radius
        if hue >= 180 :
            longitude = hue - 360
        else:
            longitude = hue
        latitude = math.degrees(math.asin(z)) 
        # radius = ((Labuv[0]-1/2)**2 + Labuv[6]**2 + Labuv[7]**2)**(1.0/2.0)
        return (latitude, longitude, radius)

    @staticmethod
    def rgb_GEOluv(rgb, profile, illuminant):
        XYZ = CVC.rgb_XYZ(rgb, profile)
        Labuv = CVC.XYZ_Labuv(XYZ, illuminant)
        GeoLuv = CVC.Labuv_GeoLuv(Labuv)
        return GeoLuv

    @staticmethod
    def rgb_GEOlab(rgb, profile, illuminant):
        XYZ = CVC.rgb_XYZ(rgb, profile)
        Labuv = CVC.XYZ_Labuv(XYZ, illuminant)
        GeoLab = CVC.Labuv_GeoLab(Labuv)
        return GeoLab
    