import ipywidgets as widgets
from IPython.display import display, clear_output
import json
from ipywidgets import Layout
layout = Layout(width='100%')
# Data dictionary
data_dict = {}

# Widgets
key_widget = widgets.Text(description='Prefix:')
label_widget = widgets.Text(description='Name')
effect_ground_pop_widget = widgets.Dropdown(options=['ref', '+', '-'], description='Effect of the mutation (if known)'
                                                                                   'in the ground state population.'
                                                                                   'ref: reference to contrast to,'
                                                                                   '+: increased population,'
                                                                                   '- decreased population.')
effect_alt1_pop_widget = widgets.Dropdown(options=['ref', '+', '-'], description='Effect of the mutation (if known)'
                                                                                   'in the most frequent alternative'
                                                                                   'state population.'
                                                                                   'ref: reference to contrast to,'
                                                                                   '+: increased population,'
                                                                                   '- decreased population.')
rank_widget = widgets.IntText(description='Rank: 0')

output = widgets.Output()

# Add to dictionary button
def add_to_dict(btn):
    with output:
        clear_output(wait=True)
        # @markdown - `sequence` Specify protein sequence to be modelled.
        key = key_widget.value
        # @markdown - `sequence` Specify protein sequence to be modelled.
        label = label_widget.value
        effect_ground_pop = effect_ground_pop_widget.value
        effect_alt1_pop = effect_alt1_pop_widget.value
        rank = rank_widget.value

        data_dict[key] = {
            "label": label,
            "effect": {"ground_pop": effect_ground_pop, "alt1_pop": effect_alt1_pop},
            "rank": rank
        }
        print(f"Added {key} to dictionary!")

add_button = widgets.Button(description="Add to Dictionary")
add_button.on_click(add_to_dict)

# Save to JSON button
def save_to_json(btn):
    with output:
        with open('data.json', 'w') as file:
            json.dump(data_dict, file)
        print("Saved to data.json!")

save_button = widgets.Button(description="Save to JSON")
save_button.on_click(save_to_json)

# Display all widgets
display(key_widget, label_widget, effect_ground_pop_widget, effect_alt1_pop_widget, rank_widget, add_button,
        save_button, output)