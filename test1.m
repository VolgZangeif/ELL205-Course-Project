[M,O] = MasterShareCreation(H2,S1,2);
retrieved = SecretImageRetrieval(H2,2,O);
%imshow(retrieved)

retrieved_comp = SecretImageRetrieval(H2_comp,2,O);
%imshow(retrieved_comp)