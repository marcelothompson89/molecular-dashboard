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
    page_title="Dashboard Molecular por País",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Función para cargar datos
@st.cache_data
def load_data():
    try:
        moleculas_compartidas = pd.read_excel('Resumen_similitud_moleculas.xlsx', 
                                            sheet_name='Moleculas compartidas')
        moleculas_unicas = pd.read_excel('Resumen_similitud_moleculas.xlsx', 
                                       sheet_name='Moleculas únicas por país')
        return moleculas_compartidas, moleculas_unicas
    except Exception as e:
        st.error(f"Error al cargar los datos: {e}")
        return None, None

# Función para filtrar datos excluyendo países
def filtrar_datos_excluidos(df_compartidas, df_unicas, paises_excluidos):
    """Filtra los dataframes excluyendo los países especificados"""
    
    # Filtrar moléculas únicas
    df_unicas_filtrado = df_unicas[~df_unicas['País'].isin(paises_excluidos)].copy()
    
    # Filtrar moléculas compartidas
    df_compartidas_filtrado = df_compartidas.copy()
    
    # Filtrar filas que solo contengan países excluidos
    indices_a_eliminar = []
    for idx, row in df_compartidas_filtrado.iterrows():
        paises_en_fila = [p.strip() for p in row['Países'].split(',')]
        paises_restantes = [p for p in paises_en_fila if p not in paises_excluidos]
        
        if len(paises_restantes) == 0:
            # Si no quedan países, eliminar la fila
            indices_a_eliminar.append(idx)
        elif len(paises_restantes) == 1:
            # Si solo queda un país, también eliminar (ya no es compartida)
            indices_a_eliminar.append(idx)
        else:
            # Actualizar la columna Países para reflejar solo los países no excluidos
            df_compartidas_filtrado.at[idx, 'Países'] = ', '.join(paises_restantes)
    
    df_compartidas_filtrado = df_compartidas_filtrado.drop(indices_a_eliminar)
    
    return df_compartidas_filtrado, df_unicas_filtrado

# Función para obtener información del país (sin cambios)
def obtener_info_pais(pais_seleccionado, df_compartidas, df_unicas):
    """Obtiene toda la información relevante del país seleccionado"""
    
    # Moléculas únicas del país
    moleculas_unicas_pais = []
    cantidad_unicas = 0
    if pais_seleccionado in df_unicas['País'].values:
        fila_pais = df_unicas[df_unicas['País'] == pais_seleccionado].iloc[0]
        cantidad_unicas = fila_pais['Cantidad de moléculas únicas']
        moleculas_unicas_pais = fila_pais['Moléculas únicas'].split(', ') if pd.notna(fila_pais['Moléculas únicas']) else []
    
    # Moléculas compartidas por el país
    moleculas_compartidas_pais = df_compartidas[
        df_compartidas['Países'].str.contains(pais_seleccionado, na=False)
    ].copy()
    
    # Colaboraciones con otros países
    colaboraciones = {}
    moleculas_por_colaboracion = {}
    
    for _, row in moleculas_compartidas_pais.iterrows():
        paises = [p.strip() for p in row['Países'].split(',')]
        otros_paises = [p for p in paises if p != pais_seleccionado]
        
        for otro_pais in otros_paises:
            if otro_pais not in colaboraciones:
                colaboraciones[otro_pais] = 0
                moleculas_por_colaboracion[otro_pais] = []
            
            colaboraciones[otro_pais] += 1
            moleculas_por_colaboracion[otro_pais].append(row['Molécula_normalizada'])
    
    # Moléculas más compartidas del país
    moleculas_mas_compartidas = []
    for _, row in moleculas_compartidas_pais.iterrows():
        num_paises = len(row['Países'].split(','))
        otros_paises = [p.strip() for p in row['Países'].split(',') if p.strip() != pais_seleccionado]
        
        moleculas_mas_compartidas.append({
            'Molécula': row['Molécula_normalizada'],
            'Número_Países': num_paises - 1,  # Excluir el país seleccionado
            'Otros_Países': ', '.join(otros_paises),
            'Países_Lista': otros_paises
        })
    
    moleculas_mas_compartidas = sorted(moleculas_mas_compartidas, 
                                     key=lambda x: x['Número_Países'], reverse=True)
    
    return {
        'moleculas_unicas': moleculas_unicas_pais,
        'cantidad_unicas': cantidad_unicas,
        'colaboraciones': colaboraciones,
        'moleculas_por_colaboracion': moleculas_por_colaboracion,
        'moleculas_mas_compartidas': moleculas_mas_compartidas,
        'total_compartidas': len(moleculas_compartidas_pais),
        'paises_colaboradores': len(colaboraciones)
    }

# Cargar datos
df_compartidas, df_unicas = load_data()

if df_compartidas is not None and df_unicas is not None:
    
    # Obtener lista de todos los países (datos originales)
    paises_compartidas = set()
    for _, row in df_compartidas.iterrows():
        paises = [p.strip() for p in row['Países'].split(',')]
        paises_compartidas.update(paises)
    
    paises_unicas = set(df_unicas['País'].tolist())
    todos_los_paises_original = sorted(list(paises_compartidas.union(paises_unicas)))
    
    # HEADER PRINCIPAL
    st.title("🧬 Análisis Molecular por País")
    st.markdown("*Explora la información molecular centrada en cada país*")
    st.markdown("---")
    
    # NUEVA SECCIÓN: FILTROS DE EXCLUSIÓN
    st.markdown("### ⚙️ Configuración de Filtros")
    
    # Crear expander para los filtros (colapsable para no ocupar mucho espacio)
    with st.expander("🚫 Excluir países del análisis", expanded=False):
        st.markdown("Selecciona los países que deseas **excluir** del análisis completo:")
        
        # Crear columnas para organizar mejor la selección
        col_excl1, col_excl2 = st.columns([3, 1])
        
        with col_excl1:
            paises_a_excluir = st.multiselect(
                "Países a excluir:",
                options=todos_los_paises_original,
                default=[],
                help="Los países seleccionados serán eliminados de todo el análisis"
            )
        
        with col_excl2:
            if paises_a_excluir:
                st.markdown("**Países excluidos:**")
                for pais in paises_a_excluir:
                    st.write(f"❌ {pais}")
            else:
                st.info("Ningún país excluido")
        
        # Botón para limpiar filtros
        if st.button("🔄 Limpiar filtros"):
            st.rerun()
    
    # Aplicar filtros de exclusión
    if paises_a_excluir:
        df_compartidas_filtrado, df_unicas_filtrado = filtrar_datos_excluidos(
            df_compartidas, df_unicas, paises_a_excluir
        )
        
        # Mostrar información sobre el filtrado
        st.info(f"📊 **Filtros aplicados:** Se han excluido {len(paises_a_excluir)} países del análisis")
    else:
        df_compartidas_filtrado = df_compartidas.copy()
        df_unicas_filtrado = df_unicas.copy()
    
    # Recalcular lista de países disponibles después del filtro
    paises_compartidas_filtrado = set()
    for _, row in df_compartidas_filtrado.iterrows():
        paises = [p.strip() for p in row['Países'].split(',')]
        paises_compartidas_filtrado.update(paises)
    
    paises_unicas_filtrado = set(df_unicas_filtrado['País'].tolist())
    todos_los_paises = sorted(list(paises_compartidas_filtrado.union(paises_unicas_filtrado)))
    
    # Verificar si hay países disponibles
    if not todos_los_paises:
        st.error("❌ No hay países disponibles después de aplicar los filtros.")
        st.stop()
    
    st.markdown("---")
    
    # SELECCIÓN DE PAÍS PRINCIPAL
    col_header1, col_header2 = st.columns([2, 1])
    
    with col_header1:
        pais_seleccionado = st.selectbox(
            "🌍 **Selecciona el país de análisis:**",
            todos_los_paises,
            index=0
        )
    
    with col_header2:
        st.markdown(f"### 📍 Analizando: **{pais_seleccionado}**")
        if paises_a_excluir:
            st.caption(f"*Datos filtrados - {len(paises_a_excluir)} países excluidos*")
    
    # Obtener información del país seleccionado (usando datos filtrados)
    info_pais = obtener_info_pais(pais_seleccionado, df_compartidas_filtrado, df_unicas_filtrado)
    
    # MÉTRICAS PRINCIPALES
    st.markdown("### 📊 Resumen del País")
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric(
            label="🔬 Moléculas Únicas",
            value=info_pais['cantidad_unicas']
        )
    
    with col2:
        st.metric(
            label="🤝 Moléculas Compartidas", 
            value=info_pais['total_compartidas']
        )
    
    with col3:
        st.metric(
            label="🌍 Países Colaboradores",
            value=info_pais['paises_colaboradores']
        )
    
    with col4:
        total_moleculas = info_pais['cantidad_unicas'] + info_pais['total_compartidas']
        st.metric(
            label="🧬 Total Moléculas",
            value=total_moleculas
        )
    
    st.markdown("---")
    
    # SECCIÓN 1: MOLÉCULAS ÚNICAS
    st.markdown("### 🔬 Moléculas Únicas del País")
    
    if info_pais['moleculas_unicas']:
        col_unique1, col_unique2 = st.columns([1, 2])
        
        with col_unique1:
            st.info(f"**{pais_seleccionado}** tiene **{len(info_pais['moleculas_unicas'])}** moléculas únicas")
            
            # Opción para mostrar/ocultar la lista
            mostrar_lista = st.checkbox("Mostrar lista completa de moléculas únicas")
        
        with col_unique2:
            if mostrar_lista:
                # Crear DataFrame para mejor visualización
                df_unicas_display = pd.DataFrame({
                    'Molécula Única': info_pais['moleculas_unicas']
                })
                df_unicas_display.index = df_unicas_display.index + 1
                st.dataframe(df_unicas_display, height=300)
    else:
        st.warning(f"**{pais_seleccionado}** no tiene moléculas únicas registradas.")
    
    st.markdown("---")
    
    # SECCIÓN 2: COLABORACIONES CON OTROS PAÍSES
    st.markdown("### 🤝 Colaboraciones Moleculares")
    
    if info_pais['colaboraciones']:
        # Top colaboradores
        top_colaboradores = dict(sorted(info_pais['colaboraciones'].items(), 
                                      key=lambda x: x[1], reverse=True)[:10])
        
        col_collab1, col_collab2 = st.columns([2, 1])
        
        with col_collab1:
            # Gráfico de barras de colaboraciones
            fig_collab = px.bar(
                x=list(top_colaboradores.values()),
                y=list(top_colaboradores.keys()),
                orientation='h',
                title=f"Top Países que Colaboran con {pais_seleccionado}",
                labels={'x': 'Moléculas Compartidas', 'y': 'País'},
                color=list(top_colaboradores.values()),
                color_continuous_scale='viridis'
            )
            fig_collab.update_layout(height=400, showlegend=False)
            st.plotly_chart(fig_collab, use_container_width=True)
        
        with col_collab2:
            st.markdown("#### 🔍 Explorar Colaboración")
            pais_colaborador = st.selectbox(
                "Selecciona un país colaborador:",
                list(top_colaboradores.keys())
            )
            
            if pais_colaborador:
                num_moleculas = info_pais['colaboraciones'][pais_colaborador]
                st.metric(f"Moléculas con {pais_colaborador}", num_moleculas)
                
                # Mostrar algunas moléculas compartidas
                moleculas_ejemplo = info_pais['moleculas_por_colaboracion'][pais_colaborador][:5]
                st.markdown("**Ejemplos de moléculas compartidas:**")
                for i, molecula in enumerate(moleculas_ejemplo, 1):
                    st.write(f"{i}. {molecula[:60]}{'...' if len(molecula) > 60 else ''}")
                
                if len(info_pais['moleculas_por_colaboracion'][pais_colaborador]) > 5:
                    st.write(f"... y {len(info_pais['moleculas_por_colaboracion'][pais_colaborador]) - 5} más")
    else:
        st.warning(f"**{pais_seleccionado}** no tiene colaboraciones moleculares registradas.")
    
    st.markdown("---")
    
    # SECCIÓN 3: MOLÉCULAS MÁS COMPARTIDAS
    st.markdown("### 🏆 Moléculas Más Compartidas del País")
    
    if info_pais['moleculas_mas_compartidas']:
        # Top moléculas más compartidas
        top_moleculas = info_pais['moleculas_mas_compartidas'][:10]
        
        col_mol1, col_mol2 = st.columns([2, 1])
        
        with col_mol1:
            # Preparar datos para el gráfico
            df_top_mol = pd.DataFrame([
                {
                    'Molécula': mol['Molécula'][:50] + "..." if len(mol['Molécula']) > 50 else mol['Molécula'],
                    'Países': mol['Número_Países']
                }
                for mol in top_moleculas
            ])
            
            fig_mol = px.bar(
                df_top_mol,
                x='Países',
                y='Molécula',
                orientation='h',
                title=f"Moléculas de {pais_seleccionado} Más Compartidas",
                labels={'x': 'Número de Países', 'y': 'Molécula'},
                color='Países',
                color_continuous_scale='plasma'
            )
            fig_mol.update_layout(height=400, showlegend=False)
            st.plotly_chart(fig_mol, use_container_width=True)
        
        with col_mol2:
            st.markdown("#### 🔬 Detalle de Molécula")
            molecula_idx = st.selectbox(
                "Selecciona una molécula:",
                range(min(10, len(top_moleculas))),
                format_func=lambda x: f"Top {x+1}: {top_moleculas[x]['Molécula'][:30]}..."
            )
            
            mol_seleccionada = top_moleculas[molecula_idx]
            st.metric("Países que la comparten", mol_seleccionada['Número_Países'])
            
            st.markdown("**Países:**")
            paises_texto = mol_seleccionada['Otros_Países']
            st.text_area("", paises_texto, height=100, disabled=True)
    else:
        st.warning(f"**{pais_seleccionado}** no tiene moléculas compartidas registradas.")
    
    st.markdown("---")
    
    # SECCIÓN 4: EXPLORADOR DE DATOS DETALLADO
    st.markdown("### 🔍 Explorador de Datos Detallado")
    
    # Filtros avanzados
    col_filter1, col_filter2, col_filter3 = st.columns(3)
    
    with col_filter1:
        tipo_molecula = st.selectbox(
            "Tipo de molécula:",
            ["Todas", "Solo Únicas", "Solo Compartidas"]
        )
    
    with col_filter2:
        min_colaboradores = st.slider(
            "Mínimo países colaboradores:",
            min_value=1,
            max_value=max([mol['Número_Países'] for mol in info_pais['moleculas_mas_compartidas']] + [1]),
            value=1
        ) if info_pais['moleculas_mas_compartidas'] else 1
    
    with col_filter3:
        buscar_texto = st.text_input("Buscar en nombres de moléculas:")
    
    # Preparar datos filtrados
    datos_filtrados = []
    
    # Agregar moléculas únicas si corresponde
    if tipo_molecula in ["Todas", "Solo Únicas"]:
        for molecula in info_pais['moleculas_unicas']:
            if not buscar_texto or buscar_texto.lower() in molecula.lower():
                datos_filtrados.append({
                    'Molécula': molecula,
                    'Tipo': 'Única',
                    'Países_Compartida': 0,
                    'Otros_Países': 'N/A'
                })
    
    # Agregar moléculas compartidas si corresponde
    if tipo_molecula in ["Todas", "Solo Compartidas"]:
        for mol_info in info_pais['moleculas_mas_compartidas']:
            if mol_info['Número_Países'] >= min_colaboradores:
                if not buscar_texto or buscar_texto.lower() in mol_info['Molécula'].lower():
                    datos_filtrados.append({
                        'Molécula': mol_info['Molécula'],
                        'Tipo': 'Compartida',
                        'Países_Compartida': mol_info['Número_Países'],
                        'Otros_Países': mol_info['Otros_Países']
                    })
    
    # Mostrar resultados
    if datos_filtrados:
        df_filtrado = pd.DataFrame(datos_filtrados)
        st.markdown(f"**Resultados encontrados:** {len(df_filtrado)} moléculas")
        st.dataframe(df_filtrado, height=400)
        
        # Opción de descarga
        csv = df_filtrado.to_csv(index=False)
        st.download_button(
            label="💾 Descargar resultados como CSV",
            data=csv,
            file_name=f"moleculas_{pais_seleccionado.replace(' ', '_')}.csv",
            mime="text/csv"
        )
    else:
        st.info("No se encontraron moléculas con los filtros aplicados.")
    
    # SECCIÓN 5: ANÁLISIS COMPARATIVO
    st.markdown("---")
    st.markdown("### ⚖️ Análisis Comparativo")
    
    # Calcular estadísticas globales para comparar
    total_paises = len(todos_los_paises)
    
    # Estadísticas del país vs promedio
    col_comp1, col_comp2, col_comp3 = st.columns(3)
    
    with col_comp1:
        # Diversidad molecular (ratio único/total)
        total_mol_pais = info_pais['cantidad_unicas'] + info_pais['total_compartidas']
        diversidad_pais = (info_pais['cantidad_unicas'] / total_mol_pais * 100) if total_mol_pais > 0 else 0
        
        st.metric(
            label="🎯 Diversidad Molecular",
            value=f"{diversidad_pais:.1f}%",
            help="Porcentaje de moléculas únicas respecto al total"
        )
    
    with col_comp2:
        # Índice de colaboración
        indice_colaboracion = (info_pais['paises_colaboradores'] / total_paises * 100) if total_paises > 0 else 0
        
        st.metric(
            label="🤝 Índice de Colaboración",
            value=f"{indice_colaboracion:.1f}%",
            help="Porcentaje de países con los que colabora"
        )
    
    with col_comp3:
        # Promedio de países por molécula compartida
        if info_pais['moleculas_mas_compartidas']:
            promedio_paises = np.mean([mol['Número_Países'] for mol in info_pais['moleculas_mas_compartidas']])
            st.metric(
                label="📊 Promedio Colaboradores/Molécula",
                value=f"{promedio_paises:.1f}",
                help="Promedio de países con los que comparte cada molécula"
            )
        else:
            st.metric(
                label="📊 Promedio Colaboradores/Molécula",
                value="0.0"
            )

else:
    st.error("No se pudieron cargar los datos. Asegúrate de que el archivo 'Resumen_similitud_moleculas.xlsx' esté en el directorio correcto.")
    st.info("El archivo debe contener las hojas: 'Moleculas compartidas' y 'Moleculas únicas por país'")

# Footer
st.markdown("---")
st.markdown("*Dashboard molecular centrado en análisis por país* 🧬")
st.markdown("*Desarrollado para explorar colaboraciones moleculares internacionales*")