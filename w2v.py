#encoding:utf-8
import numpy as np

def stopwords_removal(corpus):
    
  stop_words = set(["i", "me", "my", "myself", "we", "our", "ours", "ourselves", "you", "your", "yours", "yourself", "yourselves", "he", "him", "his", "himself", "she", "her", "hers", "herself", "it", "its", "itself", "they", "them", "their", "theirs", "themselves", "what", "which", "who", "whom", "this", "that", "these", "those", "am", "is", "are", "was", "were", "be", "been", "being", "have", "has", "had", "having", "do", "does", "did", "doing", "a", "an", "the", "and", "but", "if", "or", "because", "as", "until", "while", "of", "at", "by", "for", "with", "about", "against", "between", "into", "through", "during", "before", "after", "above", "below", "to", "from", "up", "down", "in", "out", "on", "off", "over", "under", "again", "further", "then", "once", "here", "there", "when", "where", "why", "how", "all", "any", "both", "each", "few", "more", "most", "other", "some", "such", "no", "nor", "not", "only", "own", "same", "so", "than", "too", "very", "s", "t", "can", "will", "just", "don", "should", "now"]) 

  filtered_sentence = []

  for w in corpus.split(" "): 
      if w not in stop_words: 
        filtered_sentence.append(w) 

  return " ".join(filtered_sentence)


text = """Pharmacology is the branch of pharmaceutical sciences which is concerned with the study of drug or medication action, where a drug can be broadly defined as any man-made, natural, or endogenous (from within the body) molecule which exerts a biochemical or physiological effect on the cell, tissue, organ, or organism (sometimes the word pharmacon is used as a term to encompass these endogenous and exogenous bioactive species). More specifically, it is the study of the interactions that occur between a living organism and chemicals that affect normal or abnormal biochemical function. If substances have medicinal properties, they are considered pharmaceuticals.

The field encompasses drug composition and properties, synthesis and drug design, molecular and cellular mechanisms, organ/systems mechanisms, signal transduction/cellular communication, molecular diagnostics, interactions, chemical biology, therapy, and medical applications and antipathogenic capabilities. The two main areas of pharmacology are pharmacodynamics and pharmacokinetics. Pharmacodynamics studies the effects of a drug on biological systems, and pharmacokinetics studies the effects of biological systems on a drug. In broad terms, pharmacodynamics discusses the chemicals with biological receptors, and pharmacokinetics discusses the absorption, distribution, metabolism, and excretion (ADME) of chemicals from the biological systems. Pharmacology is not synonymous with pharmacy and the two terms are frequently confused. Pharmacology, a biomedical science, deals with the research, discovery, and characterization of chemicals which show biological effects and the elucidation of cellular and organismal function in relation to these chemicals. In contrast, pharmacy, a health services profession, is concerned with the application of the principles learned from pharmacology in its clinical settings; whether it be in a dispensing or clinical care role. In either field, the primary contrast between the two is their distinctions between direct-patient care, pharmacy practice, and the science-oriented research field, driven by pharmacology.

Etymology
The word "pharmacology" is derived from Greek φάρμακον, pharmakon, "drug, poison, spell" and -λογία, -logia "study of", "knowledge of" (cf. the etymology of pharmacy). Pharmakon is related to pharmakos, the ritualistic sacrifice or exile of a human scapegoat or victim in Ancient Greek religion.

History
The origins of clinical pharmacology date back to the Middle Ages, with pharmacognosy and Avicenna's The Canon of Medicine, Peter of Spain's Commentary on Isaac, and John of St Amand's Commentary on the Antedotary of Nicholas. Early pharmacology focused on herbalism and natural substances, mainly plant extracts. Medicines were compiled in books called pharmacopoeias. Crude drugs have been used since prehistory as a preparation of substances from natural sources. However, the active ingredient of crude drugs are not purified and the substance is adulterated with other substances.

Traditional medicine varies between cultures and may be specific to a particular culture, such as in traditional Chinese, Mongolian, Tibetan and Korean medicine. However much of this has since been regarded as pseudoscience. Pharmacological substances known as entheogens may have spiritual and religious use and historical context.

In the 17th century, the English Physician Nicholas Culpeper translated and used pharmacological texts. Culpepper detailed plants and the conditions they could treat. In the 18th century, much of clinical pharmacology was established by the work of William Withering. Pharmacology as a scientific discipline did not further advance until the mid-19th century amid the great biomedical resurgence of that period. Before the second half of the nineteenth century, the remarkable potency and specificity of the actions of drugs such as morphine, quinine and digitalis were explained vaguely and with reference to extraordinary chemical powers and affinities to certain organs or tissues. The first pharmacology department was set up by Rudolf Buchheim in 1847, in recognition of the need to understand how therapeutic drugs and poisons produced their effects. Subsequently, the first pharmacology department in England was set up in 1905 at University College London.

Pharmacology developed in the 19th century as a biomedical science that applied the principles of scientific experimentation to therapeutic contexts. The advancement of research techniques propelled pharmacological research and understanding. The development of the organ bath preparation, where tissue samples are connected to recording devices, such as a myograph, and physiological responses are recorded after drug application, allowed analysis of drugs' effects on tissues. The development of the ligand binding assay in 1945 allowed quantification of the binding affinity of drugs at chemical targets. Modern pharmacologists use techniques from genetics, molecular biology, biochemistry, and other advanced tools to transform information about molecular mechanisms and targets into therapies directed against disease, defects or pathogens, and create methods for preventative care, diagnostics, and ultimately personalized medicine.

Divisions
The discipline of pharmacology can be divided into many sub disciplines each with a specific focus.

Systems of the body
Pharmacology can also focus on specific systems comprising the body. Divisions related to bodily systems study the effects of drugs in different systems of the body. These include neuropharmacology, in the central and peripheral nervous systems; immunopharmacology in the immune system. Other divisions include cardiovascular, renal and endocrine pharmacology. Psychopharmacology, is the study of the effects of drugs on the psyche, mind and behavior, such as the behavioral effects of psychoactive drugs. It incorporates approaches and techniques from neuropharmacology, animal behavior and behavioral neuroscience, and is interested in the behavioral and neurobiological mechanisms of action of psychoactive drugs. The related field of neuropsychopharmacology focuses on the effects of drugs at the overlap between the nervous system and the psyche.

Pharmacometabolomics, also known as pharmacometabonomics, is a field which stems from metabolomics, the quantification and analysis of metabolites produced by the body. It refers to the direct measurement of metabolites in an individual's bodily fluids, in order to predict or evaluate the metabolism of pharmaceutical compounds, and to better understand the pharmacokinetic profile of a drug. Pharmacometabolomics can be applied to measure metabolite levels following the administration of a drug, in order to monitor the effects of the drug on metabolic pathways. Pharmacomicrobiomics studies the effect of microbiome variations on drug disposition, action, and toxicity. Pharmacomicrobiomics is concerned with the interaction between drugs and the gut microbiome. Pharmacogenomics is the application of genomic technologies to drug discovery and further characterization of drugs related to an organism's entire genome. For pharmacology regarding individual genes, pharmacogenetics studies how genetic variation gives rise to differing responses to drugs. Pharmacoepigenetics studies the underlying epigenetic marking patterns that lead to variation in an individual's response to medical treatment.

Clinical practice and drug discovery
Pharmacology can be applied within clinical sciences. Clinical pharmacology is the basic science of pharmacology focusing on the application of pharmacological principles and methods in the medical clinic and towards patient care and outcomes. An example of this is posology, which is the study of how medicines are dosed.

Pharmacology is closely related to toxicology. Both pharmacology and toxicology are scientific disciplines that focus on understanding the properties and actions of chemicals. However, pharmacology emphasizes the therapeutic effects of chemicals, usually drugs or compounds that could become drugs, whereas toxicology is the study of chemical's adverse effects and risk assessment.

Pharmacological knowledge is used to advise pharmacotherapy in medicine and pharmacy.

Drug discovery
Drug discovery is the field of study concerned with creating new drugs. It encompases the subfields of drug design and development. Drug discovery starts with drug design, which is the inventive process of finding new drugs. In the most basic sense, this involves the design of molecules that are complementary in shape and charge to a given biomolecular target. After a lead compound has been identified through drug discovery, drug development involves bringing the drug to the market. Drug discovery is related to pharmacoeconomics, which is the sub-discipline of health economics that considers the value of drugs Pharmacoeconomics evaluates the cost and benefits of drugs in order to guide optimal healthcare resource allocation. The techniques used for the discovery, formulation, manufacturing and quality control of drugs discovery is studied by pharmaceutical engineering, a branch of engineering. Safety pharmacology specialises in detecting and investigating potential undesirable effects of drugs.

Development of medication is a vital concern to medicine, but also has strong economical and political implications. To protect the consumer and prevent abuse, many governments regulate the manufacture, sale, and administration of medication. In the United States, the main body that regulates pharmaceuticals is the Food and Drug Administration and they enforce standards set by the United States Pharmacopoeia. In the European Union, the main body that regulates pharmaceuticals is the EMA and they enforce standards set by the European Pharmacopoeia.

The metabolic stability and the reactivity of a library of candidate drug compounds have to be assessed for drug metabolism and toxicological studies. Many methods have been proposed for quantitative predictions in drug metabolism; one example of a recent computational method is SPORCalc. If the chemical structure of a medicinal compound is altered slightly, this could slightly or dramatically alter the medicinal properties of the compound depending on the level of alteration as it relates to the structural composition of the substrate or receptor site on which it exerts its medicinal effect, a concept referred to as the structural activity relationship (SAR). This means that when a useful activity has been identified, chemists will make many similar compounds called analogues, in an attempt to maximize the desired medicinal effect(s) of the compound. This development phase can take anywhere from a few years to a decade or more and is very expensive.

These new analogues need to be developed. It needs to be determined how safe the medicine is for human consumption, its stability in the human body and the best form for delivery to the desired organ system, like tablet or aerosol. After extensive testing, which can take up to 6 years, the new medicine is ready for marketing and selling.

As a result of the long time required to develop analogues and test a new medicine and the fact that of every 5000 potential new medicines typically only one will ever reach the open market, this is an expensive way of doing things, often costing over 1 billion dollars. To recoup this outlay pharmaceutical companies may do a number of things:

Carefully research the demand for their potential new product before spending an outlay of company funds.
Obtain a patent on the new medicine preventing other companies from producing that medicine for a certain allocation of time.
The inverse benefit law describes the relationship between a drugs therapeutic benefits and its marketing.

When designing drugs, the placebo effect must be considered to assess the drug's true therapeutic value.

Drug development uses techniques from medicinal chemistry to chemically design drugs. This overlaps with the biological approach of finding targets and physiological effects.

Experimentation and analysis
Theoretical pharmacology is a field of research uses techniques from computational chemistry, and molecular mechanics. Theoretical pharmacology aims at rationalizing the relation between the observed activity of a particular drug to its structural features. It aims to find relations between structure and activity. Furthermore, on the basis of the structure theoretical pharmacology aims to predict the biological activity of new drugs based on their properties and to predict new classes of drugs. Theoretical pharmacology uses pharmacometrics, which are mathematical models of biology, pharmacology, disease, and physiology used to describe and quantify interactions between drugs with pharmacology, including beneficial effects and adverse effects. Pharmacometrics can be applied to quantify drug, disease and trial information to aid efficient drug development, regulatory decisions and rational drug treatment in patients.

Experimental pharmacology involves the study of pharmacology through bioassay, to test the efficacy and potency of a drug. Systems pharmacology or network pharmacology is the application of systems biology principles in the field of pharmacology. Pharmacoinformatics relates to the broader field of bioinformatics. Ethopharmacology (not to be confused with ethnopharmacology) relates to ethology and studies drugs in the context of animal behaviours.

Wider contexts
Pharmacology can be studied in relation to wider contexts than the physiology of individuals. For example, pharmacoepidemiology is the study of the effects of drugs in large numbers of people and relates to the broader fields of epidemiology and public health. Pharmacoenvironmentology or environmental pharmacology is a field intimately linked with ecology and public health. Human health and ecology are intimately related so environmental pharmacology studies the environmental effect of drugs and pharmaceuticals and personal care products in the environment.

Drugs may also have ethnocultural importance, so ethnopharmacology studies the ethnic and cultural aspects of pharmacology.

Emerging fields
Photopharmacology is an emerging approach in medicine in which drugs are activated and deactivated with light. The energy of light is used to change for shape and chemical properties of the drug, resulting in different biological activity. This is done to ultimately achieve control when and where drugs are active in a reversible manner, to prevent side effects and pollution of drugs into the environment.

Theory of pharmacology
The study of chemicals requires intimate knowledge of the biological system affected. With the knowledge of cell biology and biochemistry increasing, the field of pharmacology has also changed substantially. It has become possible, through molecular analysis of receptors, to design chemicals that act on specific cellular signaling or metabolic pathways by affecting sites directly on cell-surface receptors (which modulate and mediate cellular signaling pathways controlling cellular function).

Chemicals can have pharmacologically relevant properties and effects. Pharmacokinetics describes the effect of the body on the chemical (e.g. half-life and volume of distribution), and pharmacodynamics describes the chemical's effect on the body (desired or toxic).

Systems, receptors and ligands
Pharmacology is typically studied with respect to particular systems, for example endogenous neurotransmitter systems. The major systems studied in pharmacology can be categorised by their ligands and include acetylcholine, adrenaline, glutamate, GABA, dopamine, histamine, serotonin, cannabinoid and opioid.

Molecular targets in pharmacology include receptors, enzymes and membrane transport proteins. Enzymes can be targeted with enzyme inhibitors. Receptors are typically categorised based on structure and function. Major receptor types studied in pharmacology include G protein coupled receptors, ligand gated ion channels and receptor tyrosine kinases.

Pharmacodynamics
Pharmacological models include the Hill equation, Cheng-Prusoff equation and Schild regression. Pharmacological theory often investigates the binding affinity of ligands to their receptors.

Medication is said to have a narrow or wide therapeutic index, certain safety factor or therapeutic window. This describes the ratio of desired effect to toxic effect. A compound with a narrow therapeutic index (close to one) exerts its desired effect at a dose close to its toxic dose. A compound with a wide therapeutic index (greater than five) exerts its desired effect at a dose substantially below its toxic dose. Those with a narrow margin are more difficult to dose and administer, and may require therapeutic drug monitoring (examples are warfarin, some antiepileptics, aminoglycoside antibiotics). Most anti-cancer drugs have a narrow therapeutic margin: toxic side-effects are almost always encountered at doses used to kill tumors.

The effect of drugs can be described with Loewe additivity.

Pharmacokinetics
Pharmacokinetics is the study of the bodily absorption, distribution, metabolism, and excretion of drugs.

When describing the pharmacokinetic properties of the chemical that is the active ingredient or active pharmaceutical ingredient (API), pharmacologists are often interested in L-ADME:

Liberation – How is the API disintegrated (for solid oral forms (breaking down into smaller particles)), dispersed, or dissolved from the medication?
Absorption – How is the API absorbed (through the skin, the intestine, the oral mucosa)?
Distribution – How does the API spread through the organism?
Metabolism – Is the API converted chemically inside the body, and into which substances. Are these active (as well)? Could they be toxic?
Excretion – How is the API excreted (through the bile, urine, breath, skin)?
Drug metabolism is assessed in pharmacokinetics and is important in drug research and prescribing.

Administration, drug policy and safety
Drug policy
In the United States, the Food and Drug Administration (FDA) is responsible for creating guidelines for the approval and use of drugs. The FDA requires that all approved drugs fulfill two requirements:

The drug must be found to be effective against the disease for which it is seeking approval (where 'effective' means only that the drug performed better than placebo or competitors in at least two trials).
The drug must meet safety criteria by being subject to animal and controlled human testing.
Gaining FDA approval usually takes several years. Testing done on animals must be extensive and must include several species to help in the evaluation of both the effectiveness and toxicity of the drug. The dosage of any drug approved for use is intended to fall within a range in which the drug produces a therapeutic effect or desired outcome.

The safety and effectiveness of prescription drugs in the U.S. is regulated by the federal Prescription Drug Marketing Act of 1987.

The Medicines and Healthcare products Regulatory Agency (MHRA) has a similar role in the UK.

Medicare Part D is a prescription drug plan in the U.S.

The Prescription Drug Marketing Act (PDMA) is an act related to drug policy.

Prescription drugs are drugs regulated by legislation.

Societies and education
Societies and administration
The International Union of Basic and Clinical Pharmacology, Federation of European Pharmacological Societies and European Association for Clinical Pharmacology and Therapeutics are an organisations representing standardisation and regulation of clinical and scientific pharmacology.

Systems for medical classification of drugs with pharmaceutical codes have been developed. These include the National Drug Code (NDC), administered by Food and Drug Administration.; Drug Identification Number (DIN), administered by Health Canada under the Food and Drugs Act; Hong Kong Drug Registration, administered by the Pharmaceutical Service of the Department of Health (Hong Kong) and National Pharmaceutical Product Index in South Africa. Hierarchical systems have also been developed, including the Anatomical Therapeutic Chemical Classification System (AT, or ATC/DDD), administered by World Health Organization; Generic Product Identifier (GPI), a hierarchical classification number published by MediSpan and SNOMED, C axis. Ingredients of drugs have been categorised by Unique Ingredient Identifier.

Education
The study of pharmacology overlaps with biomedical sciences and is study of the effects of drugs on living organisms. Pharmacological research can lead to new drug discoveries, and promote a better understanding of human physiology. Students of pharmacology must have detailed working knowledge of aspects in physiology, pathology and chemistry. Modern pharmacology is interdisciplinary and relates to biophysical and computational sciences, and analytical chemistry. Whereas a pharmacy student will eventually work in a pharmacy dispensing medications, a pharmacologist will typically work within a laboratory setting. Pharmacological research is important in academic research (medical and non-medical), private industrial positions, science writing, scientific patents and law, consultation, biotech and pharmaceutical employment, the alcohol industry, food industry, forensics/law enforcement, public health, and environmental/ecological sciences. Pharmacology is often taught to pharmacy and medicine students as part of a Medical School curriculum.""".replace(
    ")", "").replace("(", "").replace(",", "").replace(".", "")

text = stopwords_removal(text)
corpus = [[word.lower() for word in text.split(" ")]]

settings = {
    'window_size': 2,  # context window +- center word
    'n':
    50,  # dimensions of word embeddings, also refer to size of hidden layer
    'epochs': 50,  # number of training epochs
    'learning_rate': 0.01  # learning rate
}


class word2vec():
    def __init__(self):
        self.n = settings['n']
        self.lr = settings['learning_rate']
        self.epochs = settings['epochs']
        self.window = settings['window_size']

    def train(self, training_data):
        # Initialising weight matrices
        # Both s1 and s2 should be randomly initialised but for this demo, we pre-determine the arrays (getW1 and getW2)
        # getW1 - shape (9x10) and getW2 - shape (10x9)
        #self.w1 = np.array(getW1)
        #self.w2 = np.array(getW2)
        self.w1 = np.random.uniform(-1, 1, (self.v_count, self.n))
        self.w2 = np.random.uniform(-1, 1, (self.n, self.v_count))

        # Cycle through each epoch
        for i in range(self.epochs):
            # Intialise loss to 0
            
            self.loss = 0
            
            # Cycle through each training sample
            # w_t = vector for target word, w_c = vectors for context words
            for w_t, w_c in training_data:
                # Forward pass - Pass in vector for target word (w_t) to get:
                # 1. predicted y using softmax (y_pred) 2. matrix of hidden layer (h) 3. output layer before softmax (u)
                y_pred, h, u = self.forward_pass(w_t)

                EI = np.sum([np.subtract(y_pred, word) for word in w_c],
                            axis=0)

            # Backpropagation
            # We use SGD to backpropagate errors - calculate loss on the output layer
            self.backprop(EI, h, w_t)

            # Calculate loss
            # There are 2 parts to the loss function
            # Part 1: -ve sum of all the output +
            # Part 2: length of context words * log of sum for all elements (exponential-ed) in the output layer before softmax (u)
            # Note: word.index(1) returns the index in the context word vector with value 1
            # Note: u[word.index(1)] returns the value of the output layer before softmax
            self.loss += -np.sum([u[word.index(1)]
                                  for word in w_c]) + len(w_c) * np.log(
                                      np.sum(np.exp(u)))

            print("epoch: %d" % i, "loss: %d " % self.loss)
          
        print('Epoch:', i, "Loss:", self.loss)

    def backprop(self, e, h, x):
        # https://docs.scipy.org/doc/numpy-1.15.1/reference/generated/numpy.outer.html
        # Column vector EI represents row-wise sum of prediction errors across each context word for the current center word
        # Going backwards, we need to take derivative of E with respect of w2
        # h - shape 10x1, e - shape 9x1, dl_dw2 - shape 10x9
        dl_dw2 = np.outer(h, e)
        # x - shape 1x8, w2 - 5x8, e.T - 8x1
        # x - 1x8, np.dot() - 5x1, dl_dw1 - 8x5
        dl_dw1 = np.outer(x, np.dot(self.w2, e.T))
        # Update weights
        self.w1 = self.w1 - (self.lr * dl_dw1)
        self.w2 = self.w2 - (self.lr * dl_dw2)

    def word_vec(self, word):
        w_index = self.word_index[word]
        v_w = self.w1[w_index]
        
        return v_w

    def vec_sim(self, word, top_n):
      v_w1 = self.word_vec(word)
      self.vec_sim_vec(v_w1, top_n)

    def vec_sim_vec(self, v_w1, top_n):
        word_sim = {}

        for i in range(self.v_count):
            # Find the similary score for each word in vocab
            v_w2 = self.w1[i]
            theta_sum = np.dot(v_w1, v_w2)
            theta_den = np.linalg.norm(v_w1) * np.linalg.norm(v_w2)
            theta = theta_sum / theta_den

            word = self.index_word[i]
            word_sim[word] = theta

        words_sorted = sorted(
            word_sim.items(), key=lambda kv: kv[1], reverse=True)

        for word, sim in words_sorted[:top_n]:
            print(word, sim)

    def generate_training_data(self, settings, corpus):
        # Find unique word counts using dictonary
        word_counts = dict()
        for row in corpus:
            for word in row:
                try:
                    word_counts[word] += 1
                except KeyError:
                    word_counts[word] = 1
        ## How many unique words in vocab? 9
        self.v_count = len(word_counts.keys())
        # Generate Lookup Dictionaries (vocab)
        self.words_list = list(word_counts.keys())
        # Generate word:index
        self.word_index = dict(
            (word, i) for i, word in enumerate(self.words_list))
        # Generate index:word
        self.index_word = dict(
            (i, word) for i, word in enumerate(self.words_list))

        training_data = []
        # Cycle through each sentence in corpus
        for sentence in corpus:
            sent_len = len(sentence)
            # Cycle through each word in sentence
            for i, word in enumerate(sentence):
                # Convert target word to one-hot
                w_target = self.word2onehot(sentence[i])
                # Cycle through context window
                w_context = []
                # Note: window_size 2 will have range of 5 values
                for j in range(i - self.window, i + self.window + 1):
                    # Criteria for context word
                    # 1. Target word cannot be context word (j != i)
                    # 2. Index must be greater or equal than 0 (j >= 0) - if not list index out of range
                    # 3. Index must be less or equal than length of sentence (j <= sent_len-1) - if not list index out of range
                    if j != i and j <= sent_len - 1 and j >= 0:
                        # Append the one-hot representation of word to w_context
                        w_context.append(self.word2onehot(sentence[j]))
                        # print(sentence[i], sentence[j])
                        # training_data contains a one-hot representation of the target word and context words
                training_data.append([w_target, w_context])
        return np.array(training_data)


    def word2onehot(self, word):
        # word_vec - initialise a blank vector
        word_vec = [0 for i in range(0, self.v_count)
                    ]  # Alternative - np.zeros(self.v_count)
        # Get ID of word from word_index
        word_index = self.word_index[word]
        # Change value from 0 to 1 according to ID of the word
        word_vec[word_index] = 1
        return word_vec

    def forward_pass(self, x):
        # x is one-hot vector for target word, shape - 9x1
        # Run through first matrix (w1) to get hidden layer - 10x9 dot 9x1 gives us 10x1
        h = np.dot(self.w1.T, x)
        # Dot product hidden layer with second matrix (w2) - 9x10 dot 10x1 gives us 9x1
        u = np.dot(self.w2.T, h)
        # Run 1x9 through softmax to force each element to range of [0, 1] - 1x8
        y_c = self.softmax(u)
        return y_c, h, u

    def softmax(self, x):
        e_x = np.exp(x - np.max(x))
        return e_x / e_x.sum(axis=0)


def main():
  # Initialise object
  w2v = word2vec()
  # Numpy ndarray with one-hot representation for [target_word, context_words]
  training_data = w2v.generate_training_data(settings, corpus)
  w2v.train(training_data)
  w2v.vec_sim("pharmacology", 5)
  termo = input("Insira um novo termo:")
  while(termo != ""):
    w2v.vec_sim(termo, 5)
    termo = input("Insira um novo termo:")

main()
