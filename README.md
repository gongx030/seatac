# SeATAC

### [2021-01-27] Adding two branches

[development_batch_code](https://github.com/gongx030/seatac/tree/development_batch_code): In the VAE mode, each batch is represented as one hot encoding.  This is a common way of treating data from different batches in scRNA-seq analysis (for example in [this paper](https://www.biorxiv.org/content/10.1101/719856v2)).  The down-side is that this setting does not work for new data, and cannot be used for the pre-trained model.

[development_fragment_size](https://github.com/gongx030/seatac/tree/development_fragment_size): In the VAE model, the fragment size from each batch was encoded as a latent vector, concatenated with `z`, and fed into the decoder.  Since the batch effect of ATAC-seq was mostly due to the fragment size distribution, so this setting should be able to capture most of the batch effects.  Most importantly, this setting works for the pre-trained model. 

In branch [development_fragment_size](https://github.com/gongx030/seatac/tree/development_fragment_size), we will use Wald test to replace the permutation test to evaluate the significance of accessibility changes. This will significantly increase the testing speed. 
