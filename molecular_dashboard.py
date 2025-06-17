import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from collections import Counter
import re

# Configuración de la página
st.set_page_config(
    page_title="Dashboard Similitud Molecular",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Función para cargar datos
@st.cache_data
def load_data():
    try:
        # Cargar ambas hojas del Excel
        moleculas_compartidas = pd.read_excel('Resumen_similitud_moleculas.xlsx', 
                                            sheet_name='Moleculas compartidas')
        moleculas_unicas = pd.read_excel('Resumen_similitud_moleculas.xlsx', 
                                       sheet_name='Moleculas únicas por país')
        return moleculas_compartidas, moleculas_unicas
    except Exception as e:
        st.error(f"Error al cargar los datos: {e}")
        return None, None

# Función para procesar datos de países
def procesar_paises(df):
    paises_count = {}
    paises_moleculas = {}
    
    for _, row in df.iterrows():
        paises = [p.strip() for p in row['Países'].split(',')]
        molecula = row['Molécula_normalizada']
        
        for pais in paises:
            if pais not in paises_count:
                paises_count[pais] = 0
                paises_moleculas[pais] = []
            paises_count[pais] += 1
            paises_moleculas[pais].append(molecula)
    
    return paises_count, paises_moleculas

# Función para crear matriz de similitud
def crear_matriz_similitud(df):
    paises_set = set()
    for _, row in df.iterrows():
        paises = [p.strip() for p in row['Países'].split(',')]
        paises_set.update(paises)
    
    paises_list = sorted(list(paises_set))
    matriz = np.zeros((len(paises_list), len(paises_list)))
    
    for _, row in df.iterrows():
        paises = [p.strip() for p in row['Países'].split(',')]
        for i, pais1 in enumerate(paises_list):
            for j, pais2 in enumerate(paises_list):
                if pais1 in paises and pais2 in paises:
                    matriz[i][j] += 1
    
    return matriz, paises_list

# Título principal
st.title("🧬 Dashboard de Similitud Molecular entre Países")
st.markdown("---")

# Cargar datos
df_compartidas, df_unicas = load_data()

if df_compartidas is not None and df_unicas is not None:
    
    # Sidebar con información general
    st.sidebar.title("📊 Resumen General")
    st.sidebar.metric("Total Moléculas Compartidas", len(df_compartidas))
    st.sidebar.metric("Países con Moléculas Únicas", len(df_unicas))
    
    # Procesar datos
    paises_count, paises_moleculas = procesar_paises(df_compartidas)
    
    # Layout principal en columnas
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.subheader("🌍 Países con Mayor Número de Moléculas Compartidas")
        
        # Top 15 países
        top_paises = dict(sorted(paises_count.items(), key=lambda x: x[1], reverse=True)[:15])
        
        fig_bar = px.bar(
            x=list(top_paises.values()),
            y=list(top_paises.keys()),
            orientation='h',
            title="Número de Moléculas Compartidas por País",
            labels={'x': 'Número de Moléculas', 'y': 'País'},
            color=list(top_paises.values()),
            color_continuous_scale='viridis'
        )
        fig_bar.update_layout(height=600, showlegend=False)
        st.plotly_chart(fig_bar, use_container_width=True)
    
    with col2:
        st.subheader("🔬 Moléculas Únicas por País")
        
        # Gráfico de barras para moléculas únicas
        fig_unique = px.bar(
            df_unicas,
            x='Cantidad de moléculas únicas',
            y='País',
            orientation='h',
            title="Moléculas Únicas",
            color='Cantidad de moléculas únicas',
            color_continuous_scale='plasma'
        )
        fig_unique.update_layout(height=400, showlegend=False)
        st.plotly_chart(fig_unique, use_container_width=True)
        
        # Mostrar detalles de moléculas únicas
        st.subheader("📋 Detalles de Moléculas Únicas")
        pais_seleccionado = st.selectbox(
            "Selecciona un país:",
            df_unicas['País'].tolist()
        )
        
        if pais_seleccionado:
            moleculas_pais = df_unicas[df_unicas['País'] == pais_seleccionado]['Moléculas únicas'].iloc[0]
            st.text_area(
                f"Moléculas únicas de {pais_seleccionado}:",
                moleculas_pais,
                height=150
            )
    
    # Sección de matriz de similitud
    st.markdown("---")
    st.subheader("🔥 Matriz de Similitud Molecular")
    
    # Crear matriz de similitud para los top países
    top_10_paises = list(dict(sorted(paises_count.items(), key=lambda x: x[1], reverse=True)[:10]).keys())
    
    # Calcular similitud entre top países
    matriz_similitud = np.zeros((len(top_10_paises), len(top_10_paises)))
    
    for i, pais1 in enumerate(top_10_paises):
        for j, pais2 in enumerate(top_10_paises):
            if i != j:
                # Contar moléculas compartidas entre dos países específicos
                compartidas = 0
                for _, row in df_compartidas.iterrows():
                    paises_row = [p.strip() for p in row['Países'].split(',')]
                    if pais1 in paises_row and pais2 in paises_row:
                        compartidas += 1
                matriz_similitud[i][j] = compartidas
            else:
                matriz_similitud[i][j] = paises_count[pais1]
    
    # Heatmap de similitud
    fig_heatmap = px.imshow(
        matriz_similitud,
        x=top_10_paises,
        y=top_10_paises,
        color_continuous_scale='RdYlBu_r',
        title="Moléculas Compartidas entre Países (Top 10)",
        aspect="auto"
    )
    fig_heatmap.update_xaxes(side="bottom")
    fig_heatmap.update_layout(height=500)
    st.plotly_chart(fig_heatmap, use_container_width=True)
    
    # Análisis de redes
    st.markdown("---")
    st.subheader("🌐 Análisis de Colaboración Molecular")
    
    col3, col4 = st.columns(2)
    
    with col3:
        # Distribución del número de países por molécula
        paises_por_molecula = [len(row['Países'].split(',')) for _, row in df_compartidas.iterrows()]
        
        fig_dist = px.histogram(
            x=paises_por_molecula,
            nbins=20,
            title="Distribución: Número de Países por Molécula",
            labels={'x': 'Número de Países', 'y': 'Frecuencia'},
            color_discrete_sequence=['skyblue']
        )
        st.plotly_chart(fig_dist, use_container_width=True)
    
    with col4:
        # Top moléculas más compartidas
        st.subheader("🏆 Top Moléculas Más Compartidas")
        
        moleculas_mas_compartidas = []
        for _, row in df_compartidas.iterrows():
            num_paises = len(row['Países'].split(','))
            moleculas_mas_compartidas.append({
                'Molécula': row['Molécula_normalizada'][:50] + "..." if len(row['Molécula_normalizada']) > 50 else row['Molécula_normalizada'],
                'Número de Países': num_paises
            })
        
        df_top_moleculas = pd.DataFrame(moleculas_mas_compartidas)
        df_top_moleculas = df_top_moleculas.sort_values('Número de Países', ascending=False).head(10)
        
        fig_mol = px.bar(
            df_top_moleculas,
            x='Número de Países',
            y='Molécula',
            orientation='h',
            title="Top 10 Moléculas por Número de Países",
            color='Número de Países',
            color_continuous_scale='viridis'
        )
        fig_mol.update_layout(height=400, showlegend=False)
        st.plotly_chart(fig_mol, use_container_width=True)
    
    # Explorador de datos
    st.markdown("---")
    st.subheader("🔍 Explorador de Datos")
    
    # Filtros
    col5, col6 = st.columns(2)
    
    with col5:
        pais_filtro = st.selectbox(
            "Filtrar por país:",
            ['Todos'] + sorted(list(paises_count.keys()))
        )
    
    with col6:
        min_paises = st.slider(
            "Mínimo número de países por molécula:",
            min_value=2,
            max_value=max(paises_por_molecula),
            value=2
        )
    
    # Aplicar filtros
    df_filtrado = df_compartidas.copy()
    
    if pais_filtro != 'Todos':
        df_filtrado = df_filtrado[df_filtrado['Países'].str.contains(pais_filtro)]
    
    df_filtrado = df_filtrado[df_filtrado['Países'].apply(lambda x: len(x.split(',')) >= min_paises)]
    
    # Mostrar tabla filtrada
    st.subheader(f"📋 Resultados Filtrados ({len(df_filtrado)} moléculas)")
    st.dataframe(df_filtrado, height=300)
    
    # Estadísticas finales
    st.markdown("---")
    st.subheader("📈 Estadísticas Adicionales")
    
    col7, col8, col9, col10 = st.columns(4)
    
    with col7:
        st.metric("Total Países", len(paises_count))
    
    with col8:
        promedio_paises = np.mean(paises_por_molecula)
        st.metric("Promedio Países/Molécula", f"{promedio_paises:.1f}")
    
    with col9:
        max_compartidas = max(paises_por_molecula)
        st.metric("Máximo Países por Molécula", max_compartidas)
    
    with col10:
        total_unicas = df_unicas['Cantidad de moléculas únicas'].sum()
        st.metric("Total Moléculas Únicas", total_unicas)

else:
    st.error("No se pudieron cargar los datos. Asegúrate de que el archivo 'Resumen_similitud_moleculas.xlsx' esté en el directorio correcto.")
    st.info("El archivo debe contener las hojas: 'Moleculas compartidas' y 'Moleculas únicas por país'")

# Footer
st.markdown("---")
st.markdown("*Dashboard creado para análisis de similitud molecular entre países* 🧬")