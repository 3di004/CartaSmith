#Imports e Configurações ----------------------------------------------------------------------------------------------

import sympy as sp
import numpy as np
import random as rd
import matplotlib.pyplot as plt

j = 1j

#Definição da classe

class CartaSmith():
    def __init__(self, Z0, grid_color: str, figureSize: float, figdpi: int = 100):
        self.Z0 = Z0
        self.grid_color = grid_color
        self.figureSize = figureSize
        self.figdpi = figdpi

    # Método interno para criar arcos num centro C, raio r e ângulo theta.

    def ci(self, c: tuple, r: float, theta_deg: tuple = (0, 360), col: str = 'r', linewidth: float = 1, linestyle: str = '-'):
        xc, yc = c
        theta_rad: tuple = [(2*np.pi*theta_deg[i])/360 for i in range(0, 2, 1)]
        theta = np.linspace(theta_rad[0], theta_rad[1], 1000)
        x = r*np.cos(theta) + xc
        y = r*np.sin(theta) + yc
        plt.plot(x, y, c = col, ls = linestyle, lw = linewidth, zorder = 0)

    # Criador de texto ao longo de um arco

    def textCi(self, c: tuple, r: float, text: str, theta_deg: float, ajuste: float = 0, rotac: bool = True, linha: bool = False, tamL: float = 1, cor = 'default', tamDuplo:bool = False, lwd: float = 1):
        xc, yc = c
        theta_rad = (theta_deg*2*np.pi/360) - theta_deg*0.0001
        rAj = r + 0.03*r - ajuste*r
        corL = self.grid_color if cor == 'default' else cor
        x = rAj*np.cos(theta_rad) - rAj*np.cos(theta_rad)*0.02 + xc
        y = rAj*np.sin(theta_rad) + yc
        rot = (theta_deg - 90) if rotac else 0
        if not(linha):
            plt.text(x, y, text, size = 10, c = self.grid_color, rotation = rot, ha = 'center', va = 'bottom', rotation_mode = 'anchor', zorder = 0) #ha: horizontalalignment
        elif linha:
            x = r*np.cos(theta_rad*1.007)*1.005 * (1 - tamL*0.01) + xc if tamDuplo else r*np.cos(theta_rad*1.007)*1.005 + xc
            y = r*np.sin(theta_rad*1.007)*1.005 * (1 - tamL*0.01) + yc if tamDuplo else r*np.sin(theta_rad*1.007)*1.005 + yc
            xN = (x - xc) * (1 + tamL*0.02)
            yN = (y - yc) * (1 + tamL*0.02)
            plt.plot([x, xN], [y, yN], c = corL, linewidth = lwd)

    def moebius(self, X, mode: str = 'Z'): #Moebius Transformation
        if mode == 'Z':
            gamma = ((X - self.Z0)/(X + self.Z0))
        elif mode == 'Y':
            Y0 = 1/self.Z0
            gamma = ((X - Y0)/(X + Y0))
        return sp.simplify(gamma)
    
    # Função para o cálculo do coeficiente de reflexão

    def reflCoef(self, X, mode: str = 'Z') -> tuple:
        gamma = self.moebius(X)
        gammaRe = float(sp.re(gamma))
        gammaIm = float(sp.im(gamma))
        gammaC = np.sqrt(gammaRe**2 + gammaIm**2)
        ang = (np.arctan(gammaIm/gammaRe)) * (180/np.pi)
        return (gammaC, ang)
        
    # Diagrama de impedância. Retorna a figura e os axes da carta de smith

    def impedChart(self):
        fig = plt.figure(figsize = (self.figureSize, self.figureSize), dpi = self.figdpi, facecolor = "#464646")
        ax = plt.axes()

        #Círculos no y = 0 ------------------------------------------------------------------------------------------
        valorR = ['10', '5', '4', '3', '2', '1', '0,9', '0,8', '0,7', '0,6', '0,5', '0,4', '0,3', '0,2', '0,1']
        raiosR = [1/11, 1/6, 1/5, 1/4, 1/3, 1/2, 10/19,   5/9, 10/17,   5/8,   2/3,   5/7, 10/13,   5/6, 10/11]

        self.ci((0,0), 1, (0, 360), self.grid_color) #0
        for i in range(14, -1, -1):
            self.ci((1-raiosR[i],0), raiosR[i], (0, 360), self.grid_color)
        #Textos nos círculos
        ajus = [0, 0, 0, 0, 0, 0.13, 0.02, 0.015, 0.015, 0.008, 0.008, 0.01, 0, 0, 0]
        angtR = [120, 153, 148, 140, 128, 157.6, 157, 156, 154.5, 153, 151.4, 149.3, 147.1, 145, 141.5]
        for k in range(14, -1, -1):
            self.textCi((1-raiosR[k],0), raiosR[k], valorR[k], angtR[k], ajus[k])

        self.textCi((0,0), 1, '0', 180.75, -0.008)

        #Círculos no x = 1 -------------------------------------------------------------------------------------------
        valorI = ['3', '2', '1', '0,5', '0,3', '0,2', '0,1', '-0,1', '-0,2', '-0,3', '-0,5', '-1', '-2', '-3']
        valyI =  [1/3, 1/2, 1, 2, 10/3, 5, 10, -10, -5, -10/3, -2, -1, -1/2, -1/3]
        inter =  [(127, 270), (143.5, 270), (180, 270), (217, 270), (236.6, 270), (247.4, 270), (258.61, 270), (90, 101.39), (90, 112.6), (90, 123.4), (90, 143), (90, 180), (90, 216.5), (90, 233)]
        angtI =  [36.5, 53, 90, 127.7, 147.7, 158.3, 169.3, -169.3, -158.3, -147.7, -127.7, -90, -53, -36.5]

        for c in range(0, 14, 1):
            self.ci((1,valyI[c]), abs(valyI[c]), inter[c], self.grid_color)
            self.textCi((0,0), 1, valorI[c], angtI[c])

        #Círculo externo do ângulo do coeficiente de reflexão ---------------------------------------------------------

        self.ci((0,0), 1.15, (0, 360), self.grid_color)
        ang = 0
        for g in range(0, 36):
            ajuste = ang*1.006 if np.sin(ang * (np.pi/180)) < 0 else ang*1.007 if (np.sin(ang * (np.pi/180)) > 0 and np.cos(ang * (np.pi/180)) < 0) else ang
            self.textCi((0,0), 1.2, str(ang), ajuste, 0.03)
            self.textCi((0,0), 1.15, '', ang, 0, True, True)
            ang += 10

        #Retas Re e Im ------------------------------------------------------------------------------------------------
        x = np.linspace(-1, 1, 2)
        plt.plot(x, 0*x, self.grid_color, lw = 2, zorder = 0)
        plt.plot(0*x, x, self.grid_color, lw = 2, zorder = 0)

        plt.axis('off')
        plt.grid(True, ls = (0, (1, 9)))
        ax.set_xlim(-1.4, 1.4)
        ax.set_ylim(-1.4, 1.4)
        plt.title("Smith's Chart", c = 'w')
        plt.text(-1, -1.2, fr"$Z_0 = {self.Z0}\,\Omega$", c = 'w', fontsize = 11)

        plt.close()
        return fig, ax
    
    # Função para plotagem do ponto

    def plotS(self, ponto: tuple, modo:str = "Z"):
            if modo == 'Z':
                fig, ax = self.impedChart()
            fig1 = plt.figure(fig)
            ax1 = plt.axes(ax)
            ax.annotate('', xy=(1.06, 0), xycoords='axes fraction', xytext=(1.06, 1), arrowprops = dict(arrowstyle = "-", color = 'w', ls = (0, (4, 10))))
            cor = 0
            n = 1
            
            X, vec = ponto # unpacking the tuples
            gamma = self.moebius(X, modo)
            Xx = sp.re(gamma)
            Xy = sp.im(gamma)
            coefRefl = self.reflCoef(X, modo)

            cores = ['#F66B6B', '#FF61B0', '#FFDE70', '#9CDD3B', '#49ECC6', '#5A7AED', '#8E73B2', '#FC3131', '#1EFF00', '#CBC8FF']
            rd.shuffle(cores)

            plt.plot(Xx, Xy, marker = 'o', c = cores[0], zorder = 1)

            r = 1.4
            if Xx >= 0:
                theta = np.arctan(float(Xy/Xx))
            elif Xx < 0:
                theta = np.arctan(float(Xy/Xx)) + np.pi

            plt.text(r*np.cos(theta), r*np.sin(theta), fr'$\mathbf{{Z}}$', c = cores[0], fontsize = 15, fontfamily = 'Cascadia Code', zorder = 1, ha = 'center', va = 'center', rotation_mode = 'anchor')
            Yp = sp.simplify(1/X).n()
            YpRE = round(sp.re(Yp), 3)
            YpIM = round(sp.im(Yp), 3)
            plt.text(1.8, 1.3 - n/10, fr"$Z = {sp.latex(X)}\,\Omega$", c = cores[0], fontsize = 11)
            plt.text(1.8, 1.3 - (n+1)/10, fr"$Y = {sp.latex(YpRE)}{sp.latex(YpIM*j)}\,\Omega^{{-1}}$", c = cores[1], fontsize = 11)
            plt.text(1.8, 1.3 - (n+2)/10, fr"$\Gamma = {round(coefRefl[0], 3)} \angle {round(coefRefl[1], 3)}^{{\circ}}$", c = cores[2], fontsize = 11)
            Ymoeb = self.moebius(Yp, 'Y')
            plt.plot(sp.re(Ymoeb), sp.im(Ymoeb), marker = 'o', c = cores[1])
            plt.text(-r*np.cos(theta), -r*np.sin(theta), fr'$\mathbf{{Y}}$', c = cores[1], fontsize = 15, fontfamily = 'Cascadia Code', zorder = 1, ha = 'center', va = 'center', rotation_mode = 'anchor')

            #Coeficiente de onda estacionária
            self.ci((0,0), (coefRefl[0]), (0, 360), cores[3])
            s = (1 + abs(coefRefl[0]))/(1 - abs(coefRefl[0]))
            plt.plot(self.moebius(s*self.Z0), 0, marker = 'o', c = cores[3])
            plt.text(1.8, 1.3 - (n+3)/10, fr"$s = SWR = {round(s, 3)}$", c = cores[3], fontsize = 11)
            zvmax = s
            zvmin = 1/s
            plt.text(1.8, 1.3 - (n+4)/10, fr"$z_{{Vmax}} = SWR = {round(zvmax, 3)}$", c = cores[4], fontsize = 11)
            plt.text(1.8, 1.3 - (n+5)/10, fr"$z_{{Vmin}} = \dfrac{{1}}{{SWR}} = {round(zvmin, 3)}$", c = cores[5], fontsize = 11)

            self.textCi((0,0), 1.15, '', coefRefl[1], 1, True, True, 2, cores[0], True, 2.5)

            if vec:
                plt.quiver(0, 0, float(Xx), float(Xy), scale = 1, angles = 'xy', scale_units = 'xy', facecolor = cores[0], zorder = 1)
                plt.quiver(0, 0, (float(sp.re(Ymoeb))), float(sp.im(Ymoeb)), scale = 1, angles = 'xy', scale_units = 'xy', facecolor = cores[1], zorder = 1)

            fig.savefig("OutputChart.jpg", dpi = 600, bbox_inches = 'tight')
            plt.show()