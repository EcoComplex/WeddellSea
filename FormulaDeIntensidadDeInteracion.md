---
title: "Intensidad de interacción"
author: Leonardo A. Saravia
date: 18 de Nov de 2021 
output: 
  pdf_document:
  latex_engine: xelatex
editor_options: 
  chunk_output_type: console
bibliography: RegimeShifts.bib
csl: ecology.csl
---


De acuerdo a @Pawar2012 se utilizara la relación entre tamaño corporal e interacciones tróficas,  Inicialmente se propuso asumir una relación simple entre la tasa de consumo (adquisición de energía) y la tasa metabólica (uso de energía), en la cual las tasas de consumo per cápita escalan con el tamaño corporal del consumidor (m) a un exponente de aproximadamente 0,75, independientemente del taxón, el entorno o dimensionalidad [@Berlow2009;@Brown2004]. En consecuencia, las tasas de producción específicas de masa escalan como $m^{0.25}$ [@Peters1993], incluida la tasa de flujo de biomasa y las fuerzas de interacción trófica por enlace en las redes tróficas [@Otto2007; @Berlow2009; @Yodzis1992; @Brose2006; @Brose2010]. La intensidad de interacción del consumidor (predador) sobre el recurso (presa) es [@James2015; @Brose2010; @Laska1998]: 

$$q_{RC} =  \alpha x_R m_R / m_C$$ 

Donde $\alpha$ es la tasa de búsqueda, $x_R$ es la densidad media del recurso, $m_R$ es la masa corporal del recurso y $m_C$ la del consumidor. 
Ademas $\alpha / m_C$ es el coefficiente off-diagonal de las ecuaciones de Lotka-Volterra [@James2015]

La densidad media se puede estimar a partir de la biomasa corporal (Suplementary figure 2 i-j, eqn S18): 

$$x_R = X_0 m_R^{P_x}$$ 

En 2D $x_R =10^{-2.67} m_R^{-0.79}$ y en 3D $x_R=10^{-2.48} m_R^{-0.86}$

La masa corporal del recurso (Supplementary figure 2 c-d, eqn S9)

$$m_R = m_{0,R} m_C^{P_m}$$

En 2D $m_R = 10^{-2.6} m_C^{0.67}$  y en 3D $m_R=10^{-2.96} m_C ^{1.46}$


Y la tasa de búsqueda depende de la dimensionalidad del espacio de búsqueda del consumidor ( eq 3, figure 3 e-f ) asumiendo recursos escasos

si es 2D:

$$\alpha = 10^{-3.08} m_C^{0.68}$$ 

Si es 3D :

$$\alpha = 10^{-1.77} m_C ^{1.05}$$ 


Esto se puede hacer mas preciso teniendo en cuenta la estrategia de predacion que puede ser busqueda activa, filtrado (sit and wait), o pastoreo, ver material suplementario sección 1.3. aunque las diferencias entre estrategias de forrajeo son menores comparadas con las de la dimensionalidad.


## Bibliografia

