print('Hello World!')

 Add a selectbox to the sidebar:
add_selectbox = st.sidebar.selectbox(
    'How would you like to be contacted?',
    ('Email', 'Home phone', 'Mobile phone')
)

# Add a slider to the sidebar:
add_slider = st.sidebar.slider(
    'Select a range of values',
    0.0, 100.0, (25.0, 75.0)
)


button = st.sidebar.button('Compute!')
if button:
    'Starting a long computation...'

    # Add a placeholder
    latest_iteration = st.empty()
    bar = st.progress(0)

    for i in range(100):
      # Update the progress bar with each iteration.
      latest_iteration.text(f'Iteration {i+1}')
      bar.progress(i + 1)
      time.sleep(0.1)

    '...and now we\'re done!'

if os.path.exists("/app/data"):
    folders = [f.name for f in os.scandir("/app/data") if f.is_dir()]
    if folders:
        relative_path = st.sidebar.selectbox("Wähle einen Unterordner", options=folders)
        data_folder = os.path.join("/app/data", relative_path)
    else:
        st.sidebar.warning("Keine Unterordner im Datenverzeichnis gefunden.")
        data_folder = "/app/data"
else:
    st.error("Das Datenverzeichnis existiert nicht.")