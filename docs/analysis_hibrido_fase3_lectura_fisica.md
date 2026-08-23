# Curvas de carácter del híbrido: ligadura del manifold, canal 38s y lectura física de Fase 3

**Fecha**: 2026-08-22
**Autor**: Javier Aguilera
**Relevancia**: Sustituye la interpretación revocada de `analysis_hibrido_caracter_E0.md` §producción. Define la curva física correcta del híbrido Rb*(n=35)-Rb-RbCs (selección por carácter + estado de referencia sin Fermi + seguimiento adiabático), la produce a N_max=6 convergido, identifica el cruce con el canal 38s, valida el canal profundo de Fermi contra el módulo neutro, y responde con números la pregunta que motivó la Fase 3.
**Tipo**: analysis

## Resumen

La curva de ligadura del manifold n=35,l≥3 en el híbrido **no es** "el
autovalor más bajo" (eso es el umbral 37p) ni tampoco "el más bajo con
carácter punto a punto": V_Fermi^π sumerge niveles de carácter manifold
hasta cientos de GHz (−286.9 GHz en R1=800; el módulo neutro validado da
−329.2 GHz en R=800 — las figuras del paper 2016 recortan en −40 GHz).
Esos niveles son un **canal aparte**, no continuo en R1. La curva física se
define como la continuación adiabática del **autoestado más bajo con
carácter de H_A+H_mol sin Fermi** en el primer R2, seguido por solapamiento
de autovector entre puntos consecutivos. Producida a N_max=6 (convergido:
Δ(4→6)=0.013 GHz) para R1∈{600,900,1100} a₀, R2∈[500,1500] a₀:

- Las tres curvas parten del mismo objeto (peso manifold ≥99 %) y suben
  suavemente hacia el cero del manifold conforme decae el vestido molecular:
  mínimo global −23.50 / −22.51 / −23.74 GHz en R2≈500-525.
- El efecto del átomo neutro sobre la ligadura es **débil y no monótono en
  R1**: ΔE(R1,R2=500) = −0.12 / +1.07 / −0.36 GHz respecto a la referencia
  sin Fermi (~±5 % de la ligadura), oscilatorio porque atraviesa nodos
  radiales del electrón Rydberg.
- Al subir R2, la curva cruza el **umbral 38s (−20.267 GHz)** y entrega su
  carácter a ese canal: captura total desde R2≈875 para R1=600, gradual
  desde 1275 para R1=1100, sólo al final (1500) para R1=900.
- Respuesta Fase 3: el átomo neutro NO orienta ni liga al RbCs como un
  campo externo; actúa como perturbación débil (∼1 GHz) sobre la curva de
  ligadura, y como fuente de canales resonantes fuertes (niveles sumergidos
  decenas–cientos de GHz) sólo donde el electrón muestrea fuertemente el
  pseudopotencial.

## Palabras Clave

- selección por carácter (>50 % manifold)
- estado de referencia sin Fermi
- seguimiento adiabático por solapamiento
- canal 38s
- niveles sumergidos por V_Fermi
- convergencia N_max

---

## 1. Contexto

`analysis_hibrido_caracter_E0.md` revocó la lectura de las primeras curvas
E0(R2) (eran el umbral 37p) y dejó dos tareas: producir la curva por
carácter con N_max convergido, y responder la pregunta de la Fase 3
(¿el átomo neutro afecta/orienta al RbCs como lo haría un campo externo?).
Este documento cierra ambas.

## 2. Método

### 2.1 Convergencia en N_max (antes de producir)

R1=900, M_J=0, tres R2 representativos:

| N_max | dim | E(R2=500) | E(1000) | E(1500) |
|------:|----:|----------:|--------:|--------:|
| 2 | 307 | −79.47540 | −71.73432 | −71.54387 |
| 4 | 835 | −80.10430 | −71.78772 | −71.55136 |
| 6 | 1603 | −80.11733 | −71.78790 | −71.55137 |

Δ(N_max=4→6) = **0.013 GHz** < 0.1 ⇒ **N_max=6** para toda la producción.
(Estos valores corresponden al criterio puntual "más bajo con carácter"
sobre H completo; sirvieron de sonda de convergencia, no de curva final.)

### 2.2 Definición de la curva

1. **Referencia**: autoestado más bajo con peso >50 % en el manifold de
   H_A+H_mol (**sin Fermi**) evaluado en R2=500: k=49,
   E_ref = −23.3804 GHz, peso 0.9917.
2. **Arranque de cada barrido R1**: en R2=500 se toma el autoestado de H
   completo con máximo solapamiento con la referencia (ovl 0.93–0.996).
   Así las tres curvas R1 continúan EL MISMO objeto.
3. **Seguimiento**: entre puntos consecutivos de R2 (paso 25 a₀),
   j = argmax_k |⟨v_prev|v_k⟩|. No se interrumpe si pierde carácter: se
   registra el peso y se avisa (§3.3).
4. **Diagnóstico paralelo**: criterio puntual "más bajo con carácter"
   guardado siempre (K_puntual/E_puntual).

Ensamblado idéntico a `HybridNeutralPolar.hamiltonian` (verificado:
‖H−Hᵀ‖_F = 4.9e−21 Eh); V_Fermi sólo depende de R1 y se evalúa una vez por
barrido. Sanidad adicional: el método reproduce exactamente los valores de
la sonda de convergencia en R1=900.

## 3. Resultados (M_J=0, n=35, N_max=6)

### 3.1 Curvas E(R2) seguidas

| R2 (a₀) | E(R1=600) | E(R1=900) | E(R1=1100) | pesos manifold (600/900/1100) |
|--------:|----------:|----------:|-----------:|-------------------------------|
| 500 | −23.5007 | −22.3092 | −23.7400 | 0.99 / 0.99 / 0.99 |
| 600 | −20.4859 | −19.6499 | −18.7566 | 0.87 / 0.98 / 1.00 |
| 750 | −17.4569 | −17.0503 | −15.9784 | 0.99 / 1.00 / 0.98 |
| 1000 | −16.9165 | −15.6750 | −14.7291 | **0.01** / 1.00 / 0.99 |
| 1250 | −17.1548 | −15.2203 | −14.4047 | **0.01** / 1.00 / 0.76 |
| 1500 | −17.2495 | −14.3052 | −14.2820 | **0.01 / 0.06 / 0.06** |

Primer punto con peso ≤50 %: R2=875 (R1=600), 1500 (R1=900), 1275 (R1=1100).
Solapamiento mínimo entre puntos consecutivos 0.61–0.65 (los tramos bajos
de R2 cruzan bosque denso de niveles del manifold; el seguimiento aguanta).
Mínimos locales espurios tras la captura por 38s (la curva ya no es
ligadura de manifold ahí): ignorarlos.

### 3.2 Canal profundo de V_Fermi (validación cruzada con el módulo neutro)

El criterio puntual persigue estados sumergidos distintos en cada R1
(no hay objeto continuo): −23.50 (600, k=49), −286.89 (800, k=0),
−80.12 (900, k=16), −46.60 (1000, k=32), −32.05 (1100, k=43).
Aislamiento: el nivel −278.47 GHz en R1=800 existe con **Fermi solo**
(H_A+V_Fermi^π, sin H_mol), está convergido en N_max (idéntico a 4 y 6;
el rotor va de pasajero), y su composición es alta-l (l~13–31). El módulo
neutro validado (`SymmetricLinearTrimer`, paper 2016) muestra lo mismo:
suelo con pm=0.97–1.00 a −329.24/−335.91 GHz (m_l=0/1) en R=800. Conclusión:
son eigenestados reales del modelo δ-contact dentro del cloud Rydberg, un
canal propio del electrón (no artefacto del cableado del híbrido); las
figuras publicadas recortaban la ventana en −40 GHz
(`compute_trimer_curves.py`: `set_ylim(..., -40.0)`).

### 3.3 Canal 38s

El umbral 38s desnudo está a −20.267 GHz, DENTRO de la ventana de energía
de la curva [−24,−14]. Al subir R2 la curva lo cruza y le entrega el
carácter: los avisos muestran la transferencia completa
(p.ej. R1=600, R2=900: peso manifold 0.010, 38s=0.99). Por debajo de
−20.3 GHz el carácter es ≥98 % en todos los puntos: ese tramo sí es
ligadura limpia del trímero Rb*-Rb-RbCs.

## 4. Lectura física: respuesta a la pregunta de Fase 3

Pregunta: ¿puede el átomo neutro orientar o afectar al RbCs como lo haría
un campo externo?

1. **Sobre la ligadura del trímero, débilmente**: ΔE(R1,R2=500) =
   −0.12 / +1.07 / −0.36 GHz para R1=600/900/1100 frente a la referencia
   sin Fermi (ligadura molecular −23.4 GHz): efecto ≲5 %, **no monótono**
   (oscila con los nodos radiales de u_l(R1)). Un campo DC externo, en
   cambio, acopla uniformemente. El átomo no "orienta": perturba.
2. **Canales resonantes fuertes pero localizados**: donde el electrón
   muestrea fuerte el pseudopotencial aparecen niveles sumergidos
   (decenas–cientos de GHz) y ventanas de captura (38s al cruzar −20.3).
   Son espectroscopia del electrón Rydberg, no interacción RbCs-neutro
   directa.
3. **Conclusión práctica para Fase 3**: la molécula polar manda en la forma
   de la curva de ligadura (escala −14…−24 GHz, decae con R2); el átomo
   neutro introduce (i) una corrección débil ∓1 GHz oscilatoria en R1 y
   (ii) canales propios del electrón que conviene tratar aparte. No hay
   mecanismo de "orientación" análogo al campo externo en este modelo.

## 5. Limitaciones

- Paso fijo ΔR2=25 a₀: en el bosque denso cerca del cero del manifold el
  seguimiento puede saltar de rama (ovl mínimos 0.61); suficiente para
  tendencias, no para espectroscopía fina cerca de cruces.
- Los tramos con peso <50 % describen el nivel vestido 38s, no ligadura de
  manifold: no usar esos puntos como BOP.
- Pseudopotencial δ-contact: los niveles muy sumergidos (−100…−300 GHz)
  dependen sensiblemente de A_s/A_p y del corte l_max; cuantitativos sólo
  dentro de la validez del modelo compartido con el módulo neutro.
- Sólo bloque M_J=0.

## 6. Ficheros

- Script: `scripts/compute_hybrid_curves.py` (v4: referencia sin Fermi +
  tracking ininterrumpido + avisos honestos + comparación puntual)
- Datos: `plots/hybrid_neutral_polar/hybrid_caracter_R1{600,900,1100}_n35_Nmax6.npz`
  (claves: R2,E,K,Wman,Wneigh,overlap,spectrum,K_puntual,E_puntual,meta)
- Figuras: `plots/hybrid_neutral_polar/hybrid_curves_caracter_MJ0_n35.png`,
  `plots/hybrid_neutral_polar/geometria_hibrido_R1-900_R2-1000.png`
- Antecedentes: `analysis_hibrido_caracter_E0.md` (diagnóstico y revocación),
  `analysis_trimero_lineal_campo_dc.md` (módulo neutro validado)
