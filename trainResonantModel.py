import glob
import os
import pickle
import math
import gzip
import plotTools
import datetime
from timeit import default_timer as timer
from common import *
import keras

print "Imported everything in trainResonantModel.py"
inputs = [
		"jj_pt", 
        	"ll_pt",
        	"ll_M",
        	"ll_DR_l_l",
        	"jj_DR_j_j",
        	"llmetjj_DPhi_ll_jj",
        	"llmetjj_minDR_l_j",
		"llmetjj_MTformula",
		"mt2",
        	"isSF",
		"jj_M",
		"hme_h2mass_reco",
	#"hme_stddev_reco",
	#"hme_entries_reco",
	#"hme_mean_weight2_reco",
	#"hme_entry_peak_reco",
	#"hme_mostprob_offshellWmass_reco",
	#"hme_stddev_offshellWmass_reco",
	#"hme_h2mass_weight2_reco",
	#"hme_stddev_weight2_reco",
	#"hme_entries_weight2_reco"
	#"hme_mostprob_offshellWmass_weight2_reco",
	#"hme_stddev_offshellWmass_weight2_reco"
	#"hme_h2mass_reco"
        #"hme_entries_reco"
	#"hme_bin1",
        #"hme_bin2",
        #"hme_bin3",
        #"hme_bin4",
        #"hme_bin5",
        #"hme_bin6",
        #"hme_bin7",
        #"hme_bin8",
        #"hme_bin9",
        #"hme_bin10",
        #"hme_bin11",
        #"hme_bin12",
        #"hme_bin13",
        #"hme_bin14",
        #"hme_bin15",
        #"hme_bin16",
        #"hme_bin17",
        #"hme_bin18",
        #"hme_bin19",
        #"hme_bin20",
        #"hme_bin21",
        #"hme_bin22",
        #"hme_bin23",
        #"hme_bin24",
        #"hme_bin25",
        #"hme_bin26",
        #"hme_bin27",
        #"hme_bin28",
        #"hme_bin29",
        #"hme_bin30",
        #"hme_bin31",
        #"hme_bin32",
        #"hme_bin33",
        #"hme_bin34",
        #"hme_bin35",
        #"hme_bin36",
        #"hme_bin37",
        #"hme_bin38",
        #"hme_bin39",
        #"hme_bin40",
        #"hme_bin41",
        #"hme_bin42",
        #"hme_bin43",
        #"hme_bin44",
        #"hme_bin45",
        #"hme_bin46",
        #"hme_bin47",
        #"hme_bin48",
        #"hme_bin49",
        #"hme_bin50",
        #"hme_bin51",
        #"hme_bin52",
        #"hme_bin53",
        #"hme_bin54",
        #"hme_bin55",
        #"hme_bin56",
        #"hme_bin57",
        #"hme_bin58",
        #"hme_bin59",
        #"hme_bin60",
        #"hme_bin61",
        #"hme_bin62",
        #"hme_bin63",
        #"hme_bin64",
        #"hme_bin65",
        #"hme_bin66",
        #"hme_bin67",
        #"hme_bin68",
        #"hme_bin69",
        #"hme_bin70",
        #"hme_bin71",
        #"hme_bin72",
        #"hme_bin73",
        #"hme_bin74",
        #"hme_bin75",
        #"hme_bin76",
        #"hme_bin77",
        #"hme_bin78",
        #"hme_bin79",
        #"hme_bin80",
        #"hme_bin81",
        #"hme_bin82",
        #"hme_bin83",
        #"hme_bin84",
        #"hme_bin85",
        #"hme_bin86",
        #"hme_bin87",
        #"hme_bin88",
        #"hme_bin89",
        #"hme_bin90",
        #"hme_bin91",
        #"hme_bin92",
        #"hme_bin93",
        #"hme_bin94",
        #"hme_bin95",
        #"hme_bin96",
        #"hme_bin97",
        #"hme_bin98",
        #"hme_bin99"
        ]

#cut = "(91 - ll_M) > 15 && hme_h2mass_reco>=250"
cut = "(91 - ll_M) > 15"
for var in inputs:
	if var == "hme_h2mass_reco":
		cut = cut+"&& hme_h2mass_reco>=250"
		break

def lr_scheduler(epoch):
    default_lr = 0.001
    drop = 0.1
    epochs_drop = 50.0
    lr = default_lr * math.pow(drop, min(1, math.floor((1 + epoch) / epochs_drop)))
    return lr

add_mass_column = True
# add_mass_column = False
signal_masses = resonant_signal_masses

epochs = 100
suffix = str(epochs)+"epochs"
output_suffix = '{:%Y-%m-%d}_{}'.format(datetime.date.today(), suffix)
output_folder = os.path.join('hh_resonant_trained_models_kinematicwithMTandMT2MjjHME_20180206', output_suffix)

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Loading Signal and Backgrounds
print("Loading Signal and Backgrounds, cut: {}, outputfolder: {}".format(cut, output_folder))
dataset = DatasetManager(inputs, resonant_weights, cut)
dataset.load_resonant_signal(masses=signal_masses, add_mass_column=add_mass_column, fraction=1)
dataset.load_backgrounds(add_mass_column=add_mass_column)
# Split datasets into a training and testing samples (by default events are re-weighted so that S and B weight is equal). Train and test are 50%. Seed if fixed.
dataset.split()
# Those functions simple returns the tained/tested dataset/target
training_dataset, training_targets = dataset.get_training_combined_dataset_and_targets()
training_weights = dataset.get_training_combined_weights()
testing_dataset, testing_targets = dataset.get_testing_combined_dataset_and_targets()
testing_weights = dataset.get_testing_combined_weights()

output_model_filename = 'hh_resonant_trained_model.h5'
output_model_filename = os.path.join(output_folder, output_model_filename)

# callback: set of functions to be applied at given stages of the training procedure. You can use callbacks to get a view on internal states/statistics of model during training
callbacks = []
callbacks.append(keras.callbacks.ModelCheckpoint(output_model_filename, monitor='val_loss', verbose=False, save_best_only=True, mode='auto')) #Save model after every epoch
# output_logs_folder = os.path.join('hh_resonant_trained_models', 'logs', output_suffix)
# callbacks.append(keras.callbacks.TensorBoard(log_dir=output_logs_folder, histogram_freq=1, write_graph=True, write_images=False))
callbacks.append(keras.callbacks.LearningRateScheduler(lr_scheduler)) # Provide learnign rate per epoch. lr_scheduler = have to be a function of epoch.

training = True
training_time = None


# Load model
if training:
    n_inputs = len(inputs)
    if add_mass_column:
        n_inputs += 1

    # You create here the real DNN structure
    model = create_resonant_model(n_inputs)

    batch_size = 5000
    start_time = timer()
    # You do the training with the compiled model
    history = model.fit(training_dataset, training_targets, sample_weight=training_weights, batch_size=batch_size, epochs=epochs,
            verbose=True, validation_data=(testing_dataset, testing_targets, testing_weights), callbacks=callbacks)

    end_time = timer()
    training_time = datetime.timedelta(seconds=(end_time - start_time))

    save_training_parameters(output_folder, model,
            batch_size=batch_size, epochs=epochs,
            training_time=str(training_time),
            masses=signal_masses,
            with_mass_column=add_mass_column,
            inputs=inputs,
            cut=cut,
            weights=resonant_weights)

    plotTools.draw_keras_history(history, output_dir=output_folder, output_name="loss.pdf")

    # Save history
    print("Saving model training history...")
    output_history_filename = 'hh_resonant_trained_model_history.pklz'
    output_history_filename = os.path.join(output_folder, output_history_filename)
    with gzip.open(output_history_filename, 'wb') as f:
        pickle.dump(history.history, f)
    print("Done.")

model = keras.models.load_model(output_model_filename)

export_for_lwtnn(model, output_model_filename)
# Draw the inputs 
draw_resonant_training_plots(model, dataset, output_folder, split_by_mass=add_mass_column)

print("All done. Training time: %s" % str(training_time))
