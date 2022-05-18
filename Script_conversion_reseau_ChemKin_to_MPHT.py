# -*- coding: utf-8 -*-

"""

A python script that converts reaction files in the ChemKin output format, into
data files in the MPHT code input format. Don't forget to change the name of
the file you want to convert in the file parameters below before executing the
script.

"""

import copy
import os

main_dir_name = os.path.dirname(__file__)
read_uncertainties = True
add_excited_states_decay = True
add_excited_states_two_bodies_reactions = True

"""FILE PARAMETERS"""

translations_file_name = "conversion_nomenclature.txt"
reactions_file_name = "New_Curran_C2_NOX_KONNOV.inp"
uncertainties_dir_name = "Old_Chemical_Network_With_Uncertainties"
writing_dir_name = "Converted_" + reactions_file_name[:-4]

#WARNING : Information on the file structure, not meant to be changed
seps = ["<=>", "=>", "<=", "="]
keywords = ["TROE", "LOW", "SRI", "DUPLICATE", "REV", "PLOG"]
key_blocks = ["ELEMENTS", "SPECIES", "REACTIONS"]
translations_dict = {}
reverse_translations_dict = {}
# "Useless" : H2O2, cC6H6, H2, C7H8
minimum_efficiencies = ["O2", "CO", "CO2", "H2O", "CH4", "H2", "C2H6", "Ar",
                        "N2", "He", "C2H4", "cC6H6", "C7H8", "H2O2", "N2O",
                        "O3P", "NH3", "N2H4", "N2O4", "NO2", "NO", "N4S", "SO2"]

"""DECLARATION OF THE REACTION CLASS"""

class Reaction:
    """docstring"""
    def __init__(self):

        #Always defined
        self.reactants = []
        self.products = []
        self.arrhenius_coeffs = []

        #Additionnal data
        self.efficiencies = dict()
        self.translated_efficiencies = dict()
        self.stoichio = dict()
        self.troe_coeffs = []
        self.low_pressure_coeffs = []
        self.sri_coeffs = []
        self.rev_coeffs = []
        self.specific_third_bodies = []
        self.pressures = []

        #Reaction equation information
        self.is_reversible = True
        self.is_specific_third_body = False
        self.is_generic_third_body = False
        self.is_decay = False

        #Keyword information
        self.is_k_reversible = True
        self.is_troe = False
        self.is_low_pressure = False
        self.is_sri = False
        self.is_duplicate = False
        self.is_plog = False

        # Uncertainty factor on the rate coefficient ( Lognormal )
        self.F = 2
        self.F0 = 2
        # Temperature dependance of the uncertainty factor
        self.g = 100
        self.g0 = 100

"""FUNCTIONS USED"""

def verify_stoichio(species_list, species_declared):
    """docstring"""
    new_species_list = []
    for species in species_list:
        species_is_unreadable = False

        if species in species_declared or species == "M":
            new_species_list.append(species)
        else:
            try:
                stoichio = int(species[0])
            except ValueError:
                species_is_unreadable = True
                new_species_list.append(species)
            else:
                if species[1:] in species_declared or species[1:] == "M":
                    new_species_list += [species[1:]]*stoichio
                else:
                    species_is_unreadable = True
                    new_species_list.append(species)

        if species_is_unreadable:
            print("/!\ READING ERROR /!\ : THE FOLLOWING SPECIES IS UNREADABLE : " \
                  + species)
    return new_species_list

def trim_and_uncomment(file_list):
    """docstring"""
    trimmed_file = list()
    for line in file_list:
    
        #Replace tabulations with 8 whitespaces
        line = line.replace("\t", " "*8)
    
        #Remove comments
        comment_index = line.find("!")
        if comment_index != -1:
            decommented_line = line[:comment_index]
        else:
            decommented_line = line
    
        #Remove newlines and trailing and leading whitespaces
        trimmed_line = decommented_line.strip("\n ")
        if trimmed_line != "":
            trimmed_file.append(trimmed_line)
    return trimmed_file

def translate(species, reverse=False):
    """LRGP --> Olivia"""
    if reverse:
        species_dict = reverse_translations_dict
    else:
        species_dict = translations_dict
    if species in species_dict:
        return species_dict[species]
    else:
        return species

def create_species_line(species_list):
    line = ""
    for species in species_list:
        species = translate(species)
        line += " " + species + " "*(10-len(species))
    line += " "*11*(5-len(species_list))
    return line

"""----------------------------------------------------------------------------

                          READING THE INPUT FILE

----------------------------------------------------------------------------"""

"""READING THE INPUT TRANSLATION FILE"""

#Open file
with open(translations_file_name, 'r') as file:
    translations_list = file.readlines()

#Remove comments and newlines
trimmed_translations_file = trim_and_uncomment(translations_list)

#Fill the translation dictionary
for line in trimmed_translations_file:
    if "=>" in line:
       key, value = [species.strip() for species in line.split("=>")]
       translations_dict[key] = value
       reverse_translations_dict[value] = key

"""READING THE INPUT REACTIONS FILE"""

#Open file
with open(reactions_file_name, 'r') as file:
    reactions_list = file.readlines()

#Remove comments and newlines
trimmed_file = trim_and_uncomment(reactions_list)

#Identify each block
data = dict()
indexes = dict()
for key in key_blocks:
     start_index = trimmed_file.index(key)
     end_index = trimmed_file[start_index:].index("END") + start_index
     indexes[key] = [start_index, end_index]

#Go through each block
for key in key_blocks:
    data[key] = []

    #Go through the text lines corresponding to each block
    for line in trimmed_file[indexes[key][0]+1:indexes[key][1]]:

        #Read the elements and species involved
        if key != "REACTIONS":
            data[key] += line.split()

        #Read the reactions block
        else:

            #If the line corresponds to a new reaction
            if "=" in line:

                #Create a new reaction object and separate equation from data
                reaction = Reaction()
                line_data = line.split()

                #Reform the reaction equation if it had whitespaces in it
                if len(line_data) > 4:
                    line_data = ["".join(line_data[:-3])] + line_data[-3:]

                #Remove potential paranthesis from the equation
#                line_data[0] = line_data[0].replace("(", "").replace(")", "")

                #Search for separator in the equation and get involved species
                for sep in seps:
                    split_reaction = line_data[0].split(sep)

                    #Continue until correct separator found
                    if len(split_reaction) == 2:

                        #Verify if there are parenthesis around third body
                        #in each reaction term
                        for i, reaction_term in enumerate(split_reaction):
                            open_index = reaction_term.find("(+")

                            #Get rid of them
                            if open_index != -1:
                                third_body = reaction_term[open_index+2:-1]
                                split_reaction[i] = reaction_term[:open_index] \
                                + "+" + third_body

                                #Remember the specific third body if there is one
                                if third_body != "M":
                                    if third_body in data["SPECIES"]:
                                        reaction.specific_third_bodies \
                                        .append(third_body)
                                    else:
                                        print("/!\ READING ERROR /!\ : THE FOLLOWING \
                                              THIRD BODY IS UNREADABLE : " \
                                              + third_body)

                        #Get the lists of reactants and products
                        reactants, products = [reaction_term.split("+") for
                                               reaction_term in split_reaction]
                        break

                #Adjust reaction properties according to separator used
                if sep == "<=>" or sep == "=":
                    reaction.is_reversible = True
                    reaction.is_k_reversible = True
                else:
                    reaction.is_reversible = False
                    reaction.is_k_reversible = False
                    if sep == "<=":
                        reactants, products = products, reactants

                #Rewrite the stoichiometry information as repeated species
                reactants = verify_stoichio(reactants, data["SPECIES"])
                products = verify_stoichio(products, data["SPECIES"])

                #Test for generic third body in the reaction
                if "M" in reactants and "M" in products:
                    reaction.is_generic_third_body = True
                    reaction_keys = reaction.efficiencies.keys()

                    #Set default efficiencies values to 1
                    for species in data["SPECIES"]:

                        #Don't overwrite already declared efficiencies
                        if species not in reaction_keys:
                            reaction.efficiencies[species] = 1
                    reactants.remove("M")
                    products.remove("M")

                #Test for specific third body in the reaction
                third_bodies = list(set(reaction.specific_third_bodies))
                if len(third_bodies) == 1:
                    third_body = third_bodies[0]

                    #Verify if the third body appears on both sides of the
                    #reaction equation
                    nb_of_occurrences = len(reaction.specific_third_bodies)
                    if nb_of_occurrences == 2:
                        reaction.is_specific_third_body = True

                        #Set efficiencies values to 0 except for third body
                        for species in data["SPECIES"]:
                            reaction.efficiencies[species] = 0
                        reaction.efficiencies[third_body] = 1
                        reaction.specific_third_bodies = [third_body]
                        reactants.remove(third_body)
                        products.remove(third_body)

                    else:
                        print("/!\ READING ERROR /!\ : Unexpected occurrences of \
                              specific third-body " + third_body + \
                              " in reaction, detected :", nb_of_occurrences)

                elif len(third_bodies) > 1:
                    print("/!\ READING ERROR /!\ : Multiple specific third-bodies \
                          detected :", third_bodies)

                #Store information obtained from the first line
                reaction.reactants = reactants.copy()
                reaction.products = products.copy()
                reaction.arrhenius_coeffs = [float(coefficient) for coefficient \
                                             in line_data[1:]]

                data[key] += [reaction]

            #If the line corresponds to information about the previous reaction
            else:

                #Get the data between the forward slashes
                line_data = line.split("/")
                if line_data[-1] == "":
                    line_data = line_data[:-1]
                line_data = [elt.strip() for elt in line_data]

                #Analyse which keyword is used in the line
                line_is_efficiencies = True
                for keyword in keywords:
                    if keyword in line:
                        line_is_efficiencies = False
                        break

                #Get the efficiencies of the third bodies in the reaction
                if line_is_efficiencies:
                    keys, values = [], []
                    for i, elt in enumerate(line_data):
                        if i%2 == 0:
                            keys.append(elt)
                        else:
                            values.append(float(elt))
                    for coefficient_key, coefficient_value in zip(keys, values):
                        reaction.efficiencies[coefficient_key] = coefficient_value

                elif keyword == "DUPLICATE":
                    reaction.is_duplicate = True

                #Read the line according to keyword used
                else:

                    #Get the coefficients after the leading keyword
                    coeffs = [float(value) for value in line_data[-1].split()]

                    if keyword == "TROE":

                        #If the d coefficient is missing, set it to default
                        if len(coeffs) == 3:
                            coeffs += [1e18]

                        reaction.troe_coeffs = coeffs
                        reaction.is_troe = True

                    elif keyword == "LOW":

                        #If there are no troe coeffs specified, set them to 0
                        if len(reaction.troe_coeffs) == 0:
                            reaction.troe_coeffs = [0, 0, 0, 0]

                        reaction.low_pressure_coeffs = coeffs
                        reaction.is_low_pressure = True

                    elif keyword == "SRI":

                        #Unset default troe parameters
                        reaction.troe_coeffs = []

                        #If the d and e coefficient are missing,
                        #set it to 0 and 1 by default
                        if len(coeffs) == 3:
                            coeffs += [1, 0]

                        reaction.sri_coeffs = coeffs
                        reaction.is_sri = True

                    elif keyword == "REV":

                        reaction.rev_coeffs = coeffs
                        reaction.is_k_reversible = False

                    elif keyword == "PLOG":

                        #Efface les coeffs lus sur la ligne de reaction
                        if not reaction.is_plog:
                            reaction.arrhenius_coeffs = []
                        reaction.is_plog = True
                        reaction.pressures.append(coeffs[0])
                        reaction.arrhenius_coeffs.append(coeffs[1:])

"""TREAT THE READ DATA ACCORDINGLY"""

efficiency_species = []
for reaction in data["REACTIONS"]:
    #Determine efficiency species to write in the efficiencies file
    if reaction.is_generic_third_body:
        for species in reaction.efficiencies:
            if reaction.efficiencies[species] != 1 \
            and translate(species) not in efficiency_species:
                efficiency_species.append(translate(species))
    elif reaction.is_specific_third_body:
        species = reaction.specific_third_bodies[0]
        if translate(species) not in efficiency_species:
            efficiency_species.append(translate(species))

efficiency_species_to_add = [species for species in efficiency_species \
                             if species not in minimum_efficiencies]
new_efficiency_species = [species for species in minimum_efficiencies \
                          if species not in efficiency_species]
efficiency_species = minimum_efficiencies + efficiency_species_to_add

for reaction in data["REACTIONS"]:
    #Fill the translated efficiencies dictonnary
    for species in reaction.efficiencies:
        reaction.translated_efficiencies[translate(species)] = \
        reaction.efficiencies[species]

for species in new_efficiency_species:
    #Fill the new species efficiencies values that might be undefined
    for reaction in data["REACTIONS"]:
        if species not in reaction.translated_efficiencies:
            if reaction.is_generic_third_body:
                reaction.translated_efficiencies[species] = 1
            elif reaction.is_specific_third_body:
                reaction.translated_efficiencies[species] = 0

for reaction in data["REACTIONS"]:
    #Detect and tag decay reactions
    if reaction.is_generic_third_body and not reaction.is_troe \
    and len(reaction.reactants) == 1 and len(reaction.products) == 1:
        reaction.is_decay = True

for reaction in data["REACTIONS"]:
    #Convert A and C to new units ( molecules^-1 and K/J )
    Na_division_number = len(reaction.reactants)-1
    Na_rev_division_number = len(reaction.products)-1

    #If LOW is specified, then the arrhenius coeffs are high pressure
    if reaction.is_generic_third_body and not reaction.is_low_pressure:
        Na_division_number += 1
        Na_rev_division_number += 1
    if reaction.is_plog:
        for plog_arrhenius_coeffs in reaction.arrhenius_coeffs:
            plog_arrhenius_coeffs[0] /= 6.02e+23**Na_division_number
            plog_arrhenius_coeffs[2] *= 4.18/8.314
    else:
        reaction.arrhenius_coeffs[0] /= 6.02e+23**Na_division_number
        reaction.arrhenius_coeffs[2] *= 4.18/8.314

    #If REV is specified, convert the REV coeffs
    if reaction.is_reversible and not reaction.is_k_reversible:
        reaction.rev_coeffs[0] /= 6.02e+23**Na_rev_division_number
        reaction.rev_coeffs[2] *= 4.18/8.314

    #If LOW is specified, convert the LOW coeffs
    if reaction.is_low_pressure:
        reaction.low_pressure_coeffs[0] /= 6.02e+23**(len(reaction.reactants))
        reaction.low_pressure_coeffs[2] *= 4.18/8.314

reactions = []
for reaction in data["REACTIONS"]:
    #Generate the reverse reactions
    reactions.append(copy.deepcopy(reaction))
    if reaction.is_reversible and not reaction.is_k_reversible:
        reaction_copy = copy.deepcopy(reaction)
        reaction_copy.reactants, reaction_copy.products = \
        reaction_copy.products, reaction_copy.reactants
        reaction_copy.arrhenius_coeffs = reaction_copy.rev_coeffs
        reactions.append(reaction_copy)

"""READ UNCERTAINTIES"""

if read_uncertainties:
    uncertainties_dir_path = os.path.join(main_dir_name, uncertainties_dir_name)
    file_is_low_pressure = { "reaction_2_Corps_rev.dat" : False,
                             "reaction_2_Corps_irrev.dat" : False,
                             "combinaison_k0_rev.dat" : False,
                             "combinaison_k0_kinf_rev.dat" : True,
                             "combinaison_k0_kinf_irrev.dat" : True,
                             "combinaison_k0_kinf_rev_SRI.dat" : True,
                             "decompo_rev.dat" : False,
                             "decompo_irrev.dat" : False,
                             "decompo_k0_rev.dat" : False,
                             "decompo_k0_kinf_rev.dat" : True,
                             "decompo_k0_kinf_irrev.dat" : True,
                             "decompo_k0_kinf_rev_SRI.dat" : True,
                             "reaction_sansM_k0_kinf_rev.dat" : True,
                             "desexcitation_rev.dat" : False }

    file_names = list(file_is_low_pressure.keys())
    new_chem_low_reactions_infos = []
    new_chem_reactions_infos = []
    for reaction in reactions:
        if reaction.is_plog:
            new_chem_low_reactions_infos.append("PLOG FILLER")
            new_chem_reactions_infos.append("PLOG FILLER")
        else:
            reactants = [translate(species) for species in reaction.reactants]
            products = [translate(species) for species in reaction.products]
            arrhenius_coeffs = ["{:< 10.3e}".format(coeff).strip() for coeff in reaction.arrhenius_coeffs]
            low_pressure_coeffs = ["{:< 10.3e}".format(coeff).strip() for coeff in reaction.low_pressure_coeffs]
            new_chem_low_reactions_infos.append((reactants,
                                                 products,
                                                 arrhenius_coeffs,
                                                 low_pressure_coeffs))
            new_chem_reactions_infos.append((reactants,
                                             products,
                                             arrhenius_coeffs))
    for file_name in file_names:
        is_low_file = file_is_low_pressure[file_name]
        if is_low_file:
            new_chem_reactions_infos_to_replace = new_chem_low_reactions_infos
        else:
            new_chem_reactions_infos_to_replace = new_chem_reactions_infos
        with open(os.path.join(uncertainties_dir_path, file_name), 'r') as file:
            lines = file.readlines()
            for line in lines:
                reactants = line[:55].split()
                products = line[55:110].split()
                coeffs = [coeff.lower() for coeff in line[110:143].split()]
                uncertainties = [float(value) for value in line[143:165].split()]
                if is_low_file:
                    coeffs_high_p = [coeff.lower() for coeff in line[165:198].split()]
                    uncertainties_high_p = [float(value) for value in line[198:220].split()]
                    old_chem_reaction_infos = (reactants, products, coeffs_high_p, coeffs)
                else:
                    old_chem_reaction_infos = (reactants, products, coeffs)
                if uncertainties != [2, 100] or ( is_low_file and uncertainties_high_p != [2, 100] ):
                    try:
                        reaction_index = new_chem_reactions_infos_to_replace.index(old_chem_reaction_infos)
                    except ValueError:
                        continue
                    if is_low_file:
                        reactions[reaction_index].F0, reactions[reaction_index].g0 = uncertainties
                        reactions[reaction_index].F, reactions[reaction_index].g = uncertainties_high_p
                    else:
                        reactions[reaction_index].F, reactions[reaction_index].g = uncertainties

"""----------------------------------------------------------------------------

                              WRINTING THE FILES

----------------------------------------------------------------------------"""

writing_dir_path = os.path.join(main_dir_name, writing_dir_name)
if not os.path.exists(writing_dir_path):
    os.makedirs(writing_dir_path)

with open(os.path.join(writing_dir_path, "efficacites.dat"), 'w') as efficiencies_file:
    for species in efficiency_species:
        efficiencies_file.write(species + "\n")

with open(os.path.join(writing_dir_path, "reaction_2_Corps_rev.dat"), 'w')        as two_bodies_rev_file, \
     open(os.path.join(writing_dir_path, "reaction_2_Corps_rev_plog.dat"), 'w')   as two_bodies_rev_plog_file, \
     open(os.path.join(writing_dir_path, "reaction_2_Corps_irrev.dat"), 'w')      as two_bodies_irrev_file, \
     open(os.path.join(writing_dir_path, "combinaison_k0_rev.dat"), 'w')          as two_bodies_rev_coll_file, \
     open(os.path.join(writing_dir_path, "combinaison_k0_irrev.dat"), 'w')        as two_bodies_irrev_coll_file, \
     open(os.path.join(writing_dir_path, "combinaison_k0_kinf_rev.dat"), 'w')     as two_bodies_rev_troe_file, \
     open(os.path.join(writing_dir_path, "combinaison_k0_kinf_irrev.dat"), 'w')   as two_bodies_irrev_troe_file, \
     open(os.path.join(writing_dir_path, "combinaison_k0_kinf_rev_SRI.dat"), 'w') as two_bodies_rev_sri_file, \
     open(os.path.join(writing_dir_path, "decompo_rev.dat"), 'w')                 as one_body_rev_file, \
     open(os.path.join(writing_dir_path, "decompo_rev_plog.dat"), 'w')            as one_body_rev_plog_file, \
     open(os.path.join(writing_dir_path, "decompo_irrev.dat"), 'w')               as one_body_irrev_file, \
     open(os.path.join(writing_dir_path, "decompo_k0_rev.dat"), 'w')              as one_body_rev_coll_file, \
     open(os.path.join(writing_dir_path, "decompo_k0_irrev.dat"), 'w')            as one_body_irrev_coll_file, \
     open(os.path.join(writing_dir_path, "decompo_k0_kinf_rev.dat"), 'w')         as one_body_rev_troe_file, \
     open(os.path.join(writing_dir_path, "decompo_k0_kinf_irrev.dat"), 'w')       as one_body_irrev_troe_file, \
     open(os.path.join(writing_dir_path, "decompo_k0_kinf_rev_SRI.dat"), 'w')     as one_body_rev_sri_file, \
     open(os.path.join(writing_dir_path, "reaction_sansM_k0_kinf_rev.dat"), 'w')  as special_third_body_file, \
     open(os.path.join(writing_dir_path, "desexcitation_rev.dat"), 'w')           as decay_file:

    for reaction in reactions:

        #Write common info to all files
        line = create_species_line(reaction.reactants) \
             + create_species_line(reaction.products)

        if reaction.is_plog:
            pressure_list, alpha_list, beta_list, gamma_list = [], [], [], []
            for pressure, arrhenius_coeffs in zip(reaction.pressures, reaction.arrhenius_coeffs):
                alpha, beta, gamma = arrhenius_coeffs
                pressure_list.append(pressure)
                alpha_list.append(alpha)
                beta_list.append(beta)
                gamma_list.append(gamma)
            values_list = pressure_list + alpha_list + beta_list + gamma_list
            for value in values_list:
                line += " " + "{:< 10.3e}".format(value)
            line += " " + "{:< 10.3e}".format(reaction.F)
            line += " " + "{:< 10.3e}".format(reaction.g)

        else:
            #Write low pressure coeff first
            if reaction.is_low_pressure:
                for coefficient in reaction.low_pressure_coeffs:
                    line += " " + "{:< 10.3e}".format(coefficient)
                line += " " + "{:< 10.3e}".format(reaction.F0)
                line += " " + "{:< 10.3e}".format(reaction.g0)

            #Write common info to all files
            for coefficient in reaction.arrhenius_coeffs:
                line += " " + "{:< 10.3e}".format(coefficient)
            line += " " + "{:< 10.3e}".format(reaction.F)
            line += " " + "{:< 10.3e}".format(reaction.g)

            #Write everything else according to reaction type
            if reaction.is_troe or len(reaction.troe_coeffs) > 0:
                for coefficient in reaction.troe_coeffs:
                    line += " " + "{:< 10.3e}".format(coefficient)

            if reaction.is_sri:
                for coefficient in reaction.sri_coeffs:
                    line += " " + "{:< 10.3e}".format(coefficient)

            if reaction.is_specific_third_body or reaction.is_generic_third_body:
                for species in efficiency_species:
                    efficiency = reaction.translated_efficiencies[species]
                    line += " " + "{:< 10.3e}".format(efficiency)

        #Choose the correct file to write into
        if reaction.is_decay:
            file = decay_file
        elif reaction.is_specific_third_body:
            file = special_third_body_file

        #Two bodies
        elif len(reaction.reactants) == 2 or len(set(reaction.reactants)) == 2:

            #Third body
            if reaction.is_generic_third_body:

                #Troe
                if reaction.is_troe or len(reaction.troe_coeffs) > 0:
                    if reaction.is_k_reversible:
                        file = two_bodies_rev_troe_file
                    else:
                        file = two_bodies_irrev_troe_file

                #SRI
                elif reaction.is_sri:
                    file = two_bodies_rev_sri_file

                #Modified Arrhenius with third body
                else:
                    if reaction.is_k_reversible:
                        file = two_bodies_rev_coll_file
                    else:
                        file = two_bodies_irrev_coll_file

            #Simple Modified Arrhenius
            else:
                if reaction.is_k_reversible:
                    if reaction.is_plog:
                        file = two_bodies_rev_plog_file
                    else:
                        file = two_bodies_rev_file
                else:
                    file = two_bodies_irrev_file

        #One body
        elif len(reaction.reactants) == 1:

            #Third body
            if reaction.is_generic_third_body:

                #Troe
                if reaction.is_troe or len(reaction.troe_coeffs) > 0:
                    if reaction.is_k_reversible:
                        file = one_body_rev_troe_file
                    else:
                        file = one_body_irrev_troe_file

                #SRI
                elif reaction.is_sri:
                    file = one_body_rev_sri_file

                #Modified Arrhenius with third body
                else:
                    if reaction.is_k_reversible:
                        file = one_body_rev_coll_file
                    else:
                        file = one_body_irrev_coll_file

            #Simple Modified Arrhenius
            else:
                if reaction.is_k_reversible:
                    if reaction.is_plog:
                        file = one_body_rev_plog_file
                    else:
                        file = one_body_rev_file
                else:
                    file = one_body_irrev_file
        else:
             print("/!\ WRINTING ERROR /!\ : THE FOLLOWING REACTION HAS \
                   NOT ENOUGH OR TOO MANY REACTANTS : ", reaction.reactants,
                   reaction.products)

        file.write(line + "\n")

    files_to_add_dir_path = os.path.join(uncertainties_dir_path, "Files_to_add")
    if add_excited_states_decay:
        with open(os.path.join(files_to_add_dir_path, "reaction_2_Corps_rev.dat"), 'r') \
        as excited_states_two_bodies_reactions_file:
            for line in excited_states_two_bodies_reactions_file.readlines():
                two_bodies_rev_file.write(line)

    if add_excited_states_two_bodies_reactions:
        with open(os.path.join(files_to_add_dir_path, "desexcitation_rev.dat"), 'r') \
        as excited_states_decay_file:
            for line in excited_states_decay_file.readlines():
                decay_file.write(line)
