# Get all carbon sources and return the objective flux. It can be normalized by the carbon input
def get_carbon_sources(model,carbon_uptake=-10,normalize=True,original_source='EX_glc__D_e',carbon_source_list=[]):
    cobra.io.write_sbml_model(model,"tmp.xml")
    if len(carbon_source_list)==0:
        for reaction in model.reactions:
            my_id=reaction.id
            if reaction.boundary and not reaction.id.startswith("DM_"):
                for exchange_metabolite in reaction.metabolites: # Deberia haber solo un metabolito
                    if 'C' in exchange_metabolite.elements.keys():
			                n_carbons=exchange_metabolite.elements['C']
			                model_tmp=cobra.io.read_sbml_model("tmp.xml")
			                model_tmp.reactions.get_by_id(original_source).lower_bound=0 # Se modifica el uptake de la fuente original
                      if normalize:
				                model_tmp.reactions.get_by_id(my_id).lower_bound=carbon_uptake/n_carbons
                       else:
                        model_tmp.reactions.get_by_id(my_id).lower_bound=carbon_uptake
                        solution=model_tmp.optimize()
                        print("Reaction: "+reaction.id+" N carbons: "+str(n_carbons)+" Metabolite: "+exchange_metabolite.id+" Solution: "+str(solution.objective_value))
    else:
        for reaction_name in carbon_source_list:
            reaction=model.reactions.get_by_id(reaction_name)
            my_id=reaction.id
            if reaction.boundary and not reaction.id.startswith("DM_"):
                for exchange_metabolite in reaction.metabolites: # Deberia haber solo un metabolito
                    if 'C' in exchange_metabolite.elements.keys():
			                n_carbons=exchange_metabolite.elements['C']
			                model_tmp=cobra.io.read_sbml_model("tmp.xml")
			                model_tmp.reactions.get_by_id(original_source).lower_bound=0 # Se modifica el uptake de la fuente original
                      if normalize:
				                model_tmp.reactions.get_by_id(my_id).lower_bound=carbon_uptake/n_carbons
                      else:
				                model_tmp.reactions.get_by_id(my_id).lower_bound=carbon_uptake
                      solution=model_tmp.optimize()
                      print("Reaction: "+reaction.id+" N carbons: "+str(n_carbons)+" Metabolite: "+exchange_metabolite.id+" Solution: "+str(solution.objective_value))
