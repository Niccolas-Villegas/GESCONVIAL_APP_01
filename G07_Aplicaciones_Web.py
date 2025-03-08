import streamlit as st
from Lineas_influencia import *


with st.sidebar:
    st.image("Image.jpg", width=150)
    st.header('Curso: Python aplicado a la Ingeniería Estructural')
    st.subheader('Creado por: Niccolas Villegas\n')
    for _ in range(5):
        st.write('')

    r1, r2 = st.columns([1,2], gap="medium")
    with r1: "N° de tramos"
    with r2: nbeams = st.number_input("nbeams", value = 10, min_value=1, step = 1, label_visibility="collapsed")
    # NOTA: los valores deben ser del mismo formato, o es int o float

    L = []; EI = []
    for i in range(int(nbeams)):
        r1, r2, r3, r4 = st.columns([1,2,1,2], gap="small")
        with r1: st.latex(f"L_{{\\text{i+1}}}")
        with r2: Lst = st.number_input("L"+ str(i), value = 3.0, min_value=0.1, label_visibility="collapsed")
        L.append(Lst)
        with r3: st.latex(f"EI_{{\\text{i+1}}}")
        with r4: EIst = st.number_input("E"+ str(i), value = 3.0, min_value=0.1, label_visibility="collapsed")
        EI.append(EIst)

st.title('Líneas de influencia')

Type = st.radio(
        "Opciones de Análisis:",
        key="visibility",
        options=["Fuerzas Internas", "Reacciones"],
)

if Type == "Fuerzas Internas":

    Type = "Internal Forces"
    r1, r2, r3, r4  = st.columns([1,2,2,2], gap="medium")
    with r1: "Tramo N°"
    with r2: loc_sec = st.number_input("loc_sec", max_value=nbeams, min_value=1, step = 1, label_visibility="collapsed")
    if loc_sec == 1:
        dist_min = 0.1
    else:
        dist_min = 0.0
    with r3: "Distancia del inicio del tramo"
    with r4: dist_sec = st.number_input("dist_sec", max_value=L[loc_sec-1]-0.1, min_value=dist_min, step = 0.1, label_visibility="collapsed")

else:
    Type = "Reaction"
    r1, r2 = st.columns([1,2], gap="medium")
    with r1: "Apoyo N°"
    with r2: loc_sec = st.number_input("loc_sec", max_value=nbeams+1, min_value=1, step = 1, label_visibility="collapsed")
    dist_sec = 0.0


result = Influence_Lines(Type, L, EI, loc_sec, dist_sec)
fig = result.main(Type, nbeams, EI, L, loc_sec, dist_sec)

st.pyplot(fig)