Search.setIndex({"docnames": ["adder-multiplier", "chsh", "chsh_test", "dft_qpe", "hawking2", "hybrid-hhl", "intro"], "filenames": ["adder-multiplier.md", "chsh.md", "chsh_test.md", "dft_qpe.md", "hawking2.md", "hybrid-hhl.md", "intro.md"], "titles": ["Quantum adders and quantum multipliers", "The CHSH game", "The CHSH test", "Discrete Fourier Transformation and Quantum Phase Estimation", "Quantum circuit model for Blackhole evaporation", "HHL algorithm", "Topics in Quantum Computing"], "terms": {"our": [0, 2, 3, 4, 5], "first": [0, 1, 2, 3, 4, 5], "task": [0, 5], "compos": 0, "circuit": [0, 1, 2, 3, 5, 6], "perform": [0, 1, 2, 3, 5], "left": [0, 1, 2, 3, 4, 5], "right": [0, 1, 2, 3, 4, 5], "rangl": [0, 1, 2, 3, 4, 5], "mapsto": [0, 3, 5], "fix": [0, 1], "The": [0, 3, 4, 5, 6], "approach": [0, 5], "we": [0, 1, 2, 3, 4, 5], "take": [0, 1, 2, 3, 5], "us": [0, 1, 2, 3, 4, 5, 6], "discret": [0, 6], "transform": [0, 1, 4, 5, 6], "init": [0, 3, 5], "system": [0, 2, 4, 5], "state": [0, 1, 2, 3, 4, 5], "after": [0, 2, 4], "appli": [0, 1, 2, 3, 5], "dft": [0, 3, 5], "thi": [0, 1, 2, 3, 4, 5, 6], "ar": [0, 1, 2, 4, 5], "frac": [0, 1, 2, 3, 5], "1": [0, 1, 2, 3, 4, 5], "sqrt": [0, 1, 2, 3, 5], "n": [0, 3, 4, 5], "sum_": [0, 1, 3, 5], "k": [0, 1, 3, 5], "0": [0, 1, 2, 3, 4, 5], "e": [0, 1, 3, 5], "2": [0, 1, 2, 3, 4, 5], "pi": [0, 1, 2, 3, 4, 5], "iak": 0, "now": [0, 1, 2, 3, 4, 5], "consid": [0, 2, 3, 4, 5], "unitari": [0, 1, 2, 3, 4, 5], "oper": [0, 1, 2, 3, 4, 5], "u_b": 0, "i": [0, 1, 3, 4, 5, 6], "If": [0, 1, 3, 5], "get": [0, 1, 3, 4, 5], "therefor": [0, 1, 3, 4, 5], "invers": [0, 3, 4, 5], "previou": [0, 1, 3], "u_bdft": 0, "alreadi": 0, "know": [0, 1], "how": [0, 2, 5], "implement": [0, 1, 2, 5, 6], "s": [0, 1, 2, 5], "To": [0, 1, 3, 4, 5], "just": [0, 1, 3], "need": [0, 2, 3, 5], "phase": [0, 5, 6], "p": [0, 3, 5], "angl": [0, 1, 3, 4, 5], "th": [0, 3], "qubit": [0, 1, 2, 3, 4, 5, 6], "note": [0, 2], "revers": [0, 3, 5], "order": [0, 1, 3, 5, 6], "subtract": 0, "u": [0, 1, 3, 4, 5], "_b": 0, "also": [0, 1, 4, 5], "adit": 0, "modn": 0, "let": [0, 1, 2, 3, 4, 5], "write": [0, 2, 4, 6], "an": [0, 2, 3, 4, 5], "u_": [0, 3], "gate": [0, 1, 3, 4, 5], "lambda": [0, 1, 3, 4, 5], "begin": [0, 1, 3, 5], "arrai": [0, 1, 3, 4, 5], "cc": [0, 1, 3, 5], "end": [0, 1, 3, 4, 5], "lambda_": [0, 5], "otim": [0, 1, 2, 3, 5], "cdot": [0, 3, 5], "act": [0, 1, 2, 4], "k_": 0, "where": [0, 1, 3, 5], "bit": [0, 1, 3, 5, 6], "binari": [0, 3, 5], "represent": [0, 3], "given": [0, 1, 2, 3, 4, 5], "Then": [0, 1, 3], "sinc": [0, 1, 5], "pi_": 0, "final": [0, 1, 3, 4, 5], "arriv": [0, 3, 5], "from": [0, 1, 2, 3, 4, 5], "qiskit": [0, 1, 2, 4, 5, 6], "import": [0, 1, 2, 3, 4, 5], "quantumcircuit": [0, 1, 2, 3, 4, 5], "quantumregist": [0, 3, 4, 5], "classicalregist": [0, 3, 5], "quantum_info": [0, 2, 4, 5], "statevector": [0, 4, 5], "librari": [0, 5], "stateprepar": [0, 5], "numpi": [0, 1, 2, 3, 4, 5], "np": [0, 3, 4, 5], "def": [0, 1, 3, 4, 5], "dftcirc": [0, 3, 5], "make": [0, 3, 5], "rotat": [0, 1, 3, 4, 5], "rang": [0, 1, 3, 4, 5], "h": [0, 1, 2, 3, 5], "l": [0, 3, 5], "cp": [0, 3, 5], "barrier": [0, 1, 3, 4, 5], "swap": [0, 3, 4, 5], "return": [0, 1, 3, 4, 5], "below": [0, 4], "calcul": [0, 2, 3, 4, 5], "modular_sum": 0, "rais": 0, "except": [0, 4], "els": [0, 1], "qbit": 0, "cbit": 0, "modsumqc": 0, "initi": [0, 2, 3, 4], "astat": 0, "from_int": 0, "append": [0, 3, 4, 5], "measur": [0, 1, 2, 3, 5], "qc": [0, 1], "3": [0, 1, 2, 3, 4, 5], "draw": [0, 3, 4, 5], "output": [0, 3, 4, 5], "mpl": [0, 3, 4, 5], "run": [0, 1, 2, 3, 5], "simul": [0, 1, 2, 3, 4, 5], "aer": [0, 3, 5], "transpil": [0, 3, 5], "visual": [0, 3, 5], "plot_histogram": [0, 3, 5], "aer_sim": [0, 3, 5], "get_backend": [0, 3, 5], "aer_simul": [0, 3, 5], "shot": [0, 1, 2, 3, 5], "2048": [0, 3, 5], "t_modsum": 0, "result": [0, 1, 2, 3, 4, 5], "answer": [0, 1, 3, 5], "get_count": [0, 3, 5], "see": [0, 1, 3, 4, 5], "ab": [0, 3], "There": 0, "simpl": [0, 1, 3, 5], "algorithm": [0, 3, 6], "which": [0, 1, 2, 3, 4, 5], "doe": [0, 5], "time": [0, 1, 2, 4, 5], "o": 0, "simpli": [0, 5], "add": 0, "your": [0, 2], "regist": [0, 3, 5], "One": [0, 2, 3, 5], "can": [0, 1, 2, 3, 4, 5], "do": [0, 3], "obtain": [0, 1, 3, 5], "mod2": 0, "abov": [0, 1, 3], "instructionset": [0, 5], "0x7f4ba71ddde0": 0, "veri": 0, "similar": [0, 4], "section": 0, "two": [0, 1, 2, 3, 4], "gener": [0, 1, 4, 5, 6], "idea": [0, 5], "same": [0, 1, 2, 3, 4], "start": [0, 1, 2, 3, 4], "target": [0, 3, 4, 5], "matter": [0, 4], "second": [0, 1, 2, 3, 5], "understand": 0, "work": [0, 5], "must": [0, 1, 5], "condit": [0, 1, 5], "sourc": 0, "inspect": 0, "ibk": 0, "have": [0, 1, 2, 3, 4, 5], "align": [0, 1, 3, 5], "prod_": 0, "j": [0, 3, 5], "exp": [0, 3, 5], "a_": [0, 1, 2, 5], "coefici": 0, "expans": 0, "integ": [0, 3, 5], "In": [0, 1, 2, 3, 4, 5], "desir": [0, 5], "control": [0, 3, 5], "each": [0, 1, 2, 3, 4, 5], "jl": 0, "one": [0, 1, 2, 3, 4, 5, 6], "advantag": 0, "fact": [0, 1, 5], "langl": [0, 1, 2, 5], "so": [0, 1, 2, 3, 4, 5, 6], "By": [0, 1, 3, 5], "over": [0, 4], "contribuit": 0, "0l": 0, "circ": 0, "1l": 0, "all": [0, 4, 5], "modular_qsum": 0, "q1bit": 0, "q2bit": 0, "convert": [0, 2, 5], "bstate": [0, 5], "to_gat": 0, "add_regist": [0, 5], "4": [0, 1, 2, 3, 5], "5": [0, 5], "2n": 0, "auxilliari": 0, "anc": 0, "record": 0, "otimes2n": 0, "last": [0, 1, 5], "decompos": [0, 5], "r": [0, 4, 5], "b_": [0, 1, 2, 5], "i2": [0, 3], "exponenti": 0, "seen": 0, "rsl": 0, "why": 0, "ha": [0, 1, 2, 3, 5], "exactli": [0, 3, 5], "phaseg": 0, "modular_qmult": 0, "determin": [0, 3, 5], "size": [0, 5], "int": [0, 1, 2, 3, 5], "ceil": 0, "max": [0, 5], "log2": [0, 4], "abit": 0, "bbit": 0, "ancilla": 0, "multqc": 0, "c2phase": 0, "qc4": 0, "scale": [0, 3, 5], "6": 0, "0x7f4b97505660": 0, "nonloc": 1, "player": 1, "alic": [1, 2], "bob": [1, 2], "cooper": 1, "win": 1, "thei": [1, 2], "allow": [1, 5], "commun": 1, "befor": [1, 2, 4], "onc": 1, "anymor": 1, "itself": 1, "defin": [1, 3, 5], "rule": 1, "provid": [1, 5], "classic": [1, 2, 3, 5], "respond": 1, "anoth": [1, 2, 3, 4, 5], "sai": [1, 5], "x": [1, 3, 4, 5], "her": 1, "likewis": [1, 5], "y": [1, 5], "hi": 1, "b": [1, 4, 5], "follow": [1, 3, 4, 5], "satisfi": 1, "wedg": 1, "oplu": 1, "logic": 1, "conjunct": 1, "disjunct": 1, "equival": [1, 3, 5], "bitwis": 1, "addit": [1, 5], "modulo": 1, "possibl": [1, 3], "pair": [1, 2, 4], "proposit": 1, "true": [1, 2, 4, 5], "won": 1, "lost": [1, 5], "otherwis": [1, 3], "scenario": 1, "agre": [1, 3], "beforehand": 1, "whatev": 1, "receiv": [1, 2], "thu": [1, 3, 5], "incid": 1, "best": 1, "strategi": 1, "most": [1, 3], "75": 1, "success": [1, 3], "rate": 1, "what": [1, 3, 6], "interest": 1, "about": 1, "quantum": [1, 2, 5], "better": 1, "approx85": 1, "explain": [1, 3], "mean": [1, 5], "prepar": [1, 5], "bell": [1, 2], "psi": [1, 2, 3, 5], "00": [1, 2], "11": [1, 2], "For": [1, 3, 4, 5], "sake": 1, "definit": 1, "while": [1, 5], "certain": [1, 2], "hermitian": [1, 2, 5], "like": 1, "pauli": 1, "matric": [1, 2, 5], "depend": 1, "call": [1, 2, 3], "he": 1, "Their": 1, "respect": [1, 2, 4], "standard": 1, "comput": [1, 2, 3, 5], "basi": [1, 3, 5], "term": [1, 2, 5], "quantiti": 1, "choic": 1, "written": [1, 5, 6], "omega_": 1, "8": [1, 2, 3, 5], "_": [1, 2, 3, 5], "nonumb": 1, "examin": [1, 4], "expect": [1, 3, 5], "valu": [1, 2, 3, 5], "outcom": [1, 2], "other": [1, 4], "word": 1, "subindex": 1, "averag": [1, 2], "taken": 1, "13": [1, 2], "squar": 1, "ident": [1, 5], "mathbb": [1, 3, 5], "t": [1, 2, 3, 5], "properti": [1, 2], "bound": [1, 2, 5], "linear": [1, 5], "vert": [1, 5], "A": [1, 5], "leq": [1, 5], "leq4": 1, "triangl": 1, "inequ": [1, 2], "leq2": [1, 2], "leq8": 1, "rightarrow": 1, "cauchi": 1, "schwarz": 1, "go": [1, 4], "back": [1, 5], "problem": [1, 5], "hand": [1, 5], "figur": 1, "out": [1, 4], "c": [1, 2, 3, 4, 5], "matrix": [1, 3, 4, 5], "form": [1, 4], "theta": [1, 3, 5], "phi": [1, 3], "co": [1, 2, 5], "sin": [1, 2, 5], "self": 1, "adjoint": 1, "further": [1, 5], "impos": [1, 5], "2i": 1, "togeth": [1, 5], "impli": [1, 2], "odd": 1, "z": [1, 3], "abl": [1, 5], "elimin": 1, "probabilit": 1, "differ": [1, 3], "theta_": 1, "phi_": 1, "those": 1, "direct": 1, "yield": 1, "01": 1, "10": [1, 3, 5], "probabl": [1, 3, 5], "prob": 1, "reduc": [1, 3, 5], "neq": 1, "simplic": 1, "shall": 1, "assum": [1, 5], "choos": [1, 2, 5], "appropri": [1, 5], "increas": 1, "likelihood": 1, "suppos": [1, 2, 3, 5], "highest": [1, 5], "chanc": 1, "express": [1, 3, 4, 5], "9": 1, "equat": [1, 5], "alpha": [1, 4], "leftrightarrow": [1, 5], "cos2": 1, "case": [1, 3, 4], "henc": 1, "still": 1, "again": [1, 2, 3], "distinct": 1, "pm": [1, 5], "summari": 1, "propos": 1, "\u03c0": 1, "been": [1, 5], "adapt": [1, 5], "http": [1, 4], "learn": 1, "ibm": [1, 2], "com": 1, "cours": [1, 5], "basic": 1, "inform": [1, 4], "entangl": [1, 2, 4], "action": [1, 2, 4], "random": [1, 4], "randint": 1, "chsh_game": 1, "plai": 1, "arg": 1, "callabl": 1, "function": [1, 3, 4, 5], "loss": 1, "refere": 1, "randomli": 1, "decid": 1, "lose": 1, "creat": [1, 2, 3, 4, 5], "base": [1, 3, 5], "present": 1, "chsh_circuit": 1, "when": [1, 3, 4, 5], "paramet": [1, 5], "cx": [1, 2, 4, 5], "qiskit_a": 1, "primit": [1, 2], "sampler": 1, "quantum_strategi": 1, "carri": 1, "statist": [1, 5], "quasi_dist": 1, "binary_prob": 1, "list": [1, 4, 5], "kei": [1, 5], "num_gam": 1, "1000": 1, "total_scor": 1, "print": [1, 2, 3, 4, 5], "fraction": 1, "83": 1, "theorem": [2, 5], "assert": 2, "mechan": [2, 4], "incompat": 2, "local": 2, "hidden": 2, "variabl": 2, "theori": 2, "constraint": 2, "give": [2, 3], "rise": 2, "name": [2, 3, 4, 5], "author": 2, "clauser": 2, "horn": 2, "shimoni": 2, "holt": 2, "experiment": 2, "prove": 2, "particl": 2, "moreov": 2, "suffici": 2, "far": 2, "apart": 2, "cannot": 2, "influenc": 2, "hypothesi": 2, "These": 2, "correspond": [2, 3, 4, 5], "forth": 2, "alwai": 2, "pm1": [2, 5], "combin": [2, 5], "equal": 2, "instanc": [2, 3, 5], "zero": [2, 5], "pm2": 2, "happen": [2, 4], "opposit": 2, "becaus": [2, 5], "experi": 2, "describ": 2, "repeat": [2, 3, 5], "mani": 2, "12": 2, "less": [2, 3], "than": [2, 6], "sum": [2, 4, 5], "label": 2, "eq": 2, "post": 2, "game": [2, 6], "shown": 2, "violat": [2, 4], "briefli": 2, "main": 2, "argument": 2, "belong": [2, 5], "And": 2, "analog": 2, "show": [2, 4], "exce": 2, "upper": [2, 5], "deduc": 2, "real": [2, 5], "estim": [2, 5, 6], "extens": [2, 5], "unitaryg": 2, "runtim": 2, "qiskit_ibm_runtim": 2, "qiskitruntimeservic": 2, "session": 2, "option": 2, "down": 2, "a0_matrix": 2, "a1_matrix": 2, "b0_matrix": 2, "b1_matrix": 2, "class": 2, "tensor": [2, 4], "product": [2, 4], "a0": 2, "a1": 2, "b0": 2, "b1": 2, "a0b0": 2, "a0b1": 2, "a1b0": 2, "a1b1": 2, "op": 2, "next": [2, 3, 4, 5], "bellqc": 2, "number": [2, 3, 5], "backend": [2, 5], "servic": 2, "channel": 2, "ibm_quantum": 2, "q": 2, "research": 2, "fed": 2, "uni": 2, "ufabc": 2, "token": 2, "b991d7c5d1d6efcf2ec598451d82108070c330d72f3b033846e4aaecff48ae1c65575b7463772cf57e1fc8591310f30e299ed962d96f4e6de74e33d44f4db836": 2, "set": [2, 5], "overwritten": 2, "job": 2, "level": 2, "optimization_level": 2, "resilience_level": 2, "select": 2, "fewest": 2, "queue": 2, "least_busi": 2, "fals": 2, "ibm_lago": 2, "observ": 2, "1e4": 2, "f": [2, 3, 4, 5], "id": [2, 5], "job_id": 2, "statu": 2, "home": 2, "fresneda": 2, "lib": 2, "python3": 2, "site": 2, "packag": 2, "qiskit_runtime_servic": 2, "py": 2, "994": 2, "userwarn": 2, "current": 2, "dedic": 2, "warn": 2, "cmg6kyqdyqh0008rscmg": 2, "jobstatu": 2, "queu": 2, "retriev": 2, "cjmbo7cvcjlre5d4oarg": 2, "compar": [2, 4, 5], "upperbound": 2, "35101117022943606": 2, "chsh2": 2, "3767145006894097": 2, "here": 3, "ijk": 3, "j_": 3, "j_i": 3, "foral": 3, "wai": [3, 5], "i0": 3, "rational": 3, "hadamard": [3, 5], "produc": 3, "j_0": 3, "u_2": 3, "u_3": 3, "u_4": 3, "u_n": 3, "bring": [3, 5], "applic": 3, "procedur": 3, "remain": 3, "defint": 3, "qft_rotat": 3, "exit": 3, "empti": 3, "index": 3, "signific": 3, "smaller": 3, "singl": 3, "At": [3, 4], "earlier": 3, "swap_regist": 3, "qft": 3, "eigenvector": [3, 5], "eigenvalu": [3, 4, 5], "varphi": 3, "aim": 3, "find": [3, 5, 6], "necessari": [3, 5], "store": [3, 5], "_2": 3, "accur": 3, "least": 3, "epsilon": 3, "lceil": 3, "log": 3, "rceil": 3, "count": [3, 5], "relat": [3, 5], "ik": 3, "ignor": 3, "dagger": [3, 5], "tild": [3, 5], "accuraci": [3, 5], "varphi_": 3, "varphi_0": 3, "stage": 3, "would": [3, 5, 6], "made": [3, 5], "explicitli": 3, "aproxim": [3, 5], "peak": 3, "chose": 3, "eigenst": [3, 5], "exact": 3, "approxim": [3, 4, 5], "cqbit": [3, 5], "mbit": [3, 5], "repetit": [3, 5], "counting_qubit": [3, 5], "t_qpe": [3, 5], "001": 3, "times2": [3, 5], "turn": 3, "chang": 3, "distribut": [3, 4, 5], "It": [3, 5], "difficult": 3, "read": 3, "pictur": 3, "wa": 3, "frequent": 3, "most_frequ": 3, "3330078125": 3, "percent": 3, "round": [3, 5], "100": 3, "off": 3, "causal": 4, "arxiv": 4, "org": 4, "pdf": 4, "2104": 4, "14901": 4, "four": 4, "being": 4, "m": [4, 5], "arbitrari": 4, "mu": 4, "repres": 4, "fall": 4, "through": 4, "event": 4, "horizon": 4, "beta": 4, "coupl": 4, "radiat": 4, "g": 4, "pass": 4, "cnot": 4, "rangle_m": 4, "rangle_g": 4, "hawk": 4, "r_1": 4, "insid": 4, "bh": 4, "r_2": 4, "escap": 4, "process": 4, "u_h": 4, "kappa": 4, "rangle_": 4, "outsid": 4, "more": [4, 6], "physic": 4, "blackol": 4, "transfer": 4, "whole": 4, "never": 4, "thank": 4, "among": 4, "howev": [4, 5], "complet": [4, 5], "outgo": 4, "densitymatrix": [4, 5], "partial_trac": [4, 5], "random_statevector": 4, "hwqc": 4, "infal": 4, "initial_st": 4, "data": 4, "d0": 4, "d1": 4, "r0": 4, "r1": 4, "d2": 4, "between": [4, 5], "mode": 4, "d3": 4, "d4": 4, "trace": 4, "inner": 4, "p0": 4, "p1": 4, "p2": 4, "p3": 4, "p4": 4, "both": 4, "matter0": 4, "hw1": 4, "page": 4, "profil": 4, "step": [4, 5], "save": 4, "densiti": 4, "entropi": 4, "plot": 4, "curv": 4, "von": 4, "neumann": 4, "neumann_entropi": 4, "lambda_k": 4, "ln": 4, "lamda_k": 4, "exclud": 4, "eigen": 4, "linalg": [4, 5], "eigvalsh": 4, "000001": 4, "matplotlib": 4, "pyplot": 4, "plt": 4, "bo": 4, "ylabel": 4, "s_r": 4, "xlabel": 4, "suptitl": 4, "try": 4, "construct": 4, "incom": 4, "bh_evapor": 4, "genhwqc": 4, "non": [4, 5], "degre": 4, "freedom": 4, "dens1": 4, "entropies1": 4, "dens2": 4, "entropies2": 4, "check": [4, 5], "statement": 4, "dens3": 4, "entropies3": 4, "discuss": 5, "solv": 5, "ax": 5, "mathcal": 5, "complex": 5, "without": 5, "usual": [5, 6], "norm": 5, "conjug": 5, "transpos": 5, "assumpt": 5, "rescal": 5, "solut": 5, "realiz": 5, "coeffici": 5, "invert": 5, "uniqu": 5, "spectral": 5, "lie": 5, "compact": 5, "min": 5, "subset": 5, "orthonorm": 5, "v_": 5, "orthogon": 5, "project": 5, "diagon": 5, "rather": 5, "its": 5, "extend": 5, "becom": 5, "m_": 5, "altern": 5, "necess": 5, "interv": 5, "accomplish": 5, "shift": 5, "spectrum": 5, "prior": 5, "knowledg": 5, "point": 5, "achiev": 5, "radiu": 5, "rho": 5, "ani": 5, "conveni": 5, "maximum": 5, "absolut": 5, "row": 5, "equiv": 5, "infti": 5, "max_": 5, "ij": 5, "spec": 5, "omega": 5, "root": 5, "characterist": 5, "polynomi": 5, "det": 5, "tell": 5, "notic": 5, "v": 5, "replac": 5, "accord": 5, "deatlt": 5, "whose": 5, "onli": 5, "instead": 5, "ia": 5, "qpe": [5, 6], "overset": 5, "should": 5, "reciproc": 5, "appear": 5, "introduc": 5, "auxiliar": 5, "r_": 5, "arcsin": 5, "could": 5, "arcco": 5, "up": 5, "factor": 5, "fit": 5, "domain": 5, "effici": 5, "seem": 5, "ref1": 5, "done": 5, "newton": 5, "method": 5, "bissect": 5, "ref2": 5, "ref3": 5, "piecewis": 5, "adopt": [5, 6], "third": 5, "involv": 5, "multipl": 5, "ry": 5, "put": 5, "everyth": 5, "respons": 5, "known": 5, "new": 5, "shiftop": 5, "shape": 5, "inf": 5, "shifta": 5, "eigenv": 5, "eigashift": 5, "ashift": 5, "lamb": 5, "As": 5, "exampl": 5, "test": [5, 6], "dtype": 5, "complex128": 5, "shifth": 5, "eighshift": 5, "eigval": 5, "subroutin": 5, "hamiltoniang": 5, "shiftopg": 5, "rygat": 5, "controlledg": 5, "precison": 5, "qpeqc": 5, "0x7fcf401642e0": 5, "sorted_list": 5, "sort": 5, "item": 5, "phi1": 5, "phi2": 5, "lambda1": 5, "lambda2": 5, "phi1b": 5, "phi2b": 5, "9998931827716024": 5, "9998931827716022": 5, "three": 5, "encod": 5, "targetbit": 5, "hhlqc": 5, "minimum": 5, "ry1": 5, "num_ctrl_qubit": 5, "ctrl_state": 5, "ry2": 5, "0x7fcf397d43a0": 5, "analyz": 5, "custom": 5, "basis_g": 5, "phi1gat": 5, "t_hhl": 5, "histogram": 5, "close": 5, "collaps": 5, "fidel": 5, "vector": 5, "statevectorsimul": 5, "execut": 5, "state_fidel": 5, "statevector_simul": 5, "get_statevector": 5, "traced_sv": 5, "cstate": 5, "9994": 5, "excel": 5, "agreement": 5, "collect": 6, "fill": 6, "littl": 6, "mathemat": 6, "detail": 6, "onlin": 6, "materi": 6, "lack": 6, "some": 6, "clariti": 6, "convent": 6, "endian": 6, "oppos": 6, "literatur": 6, "easier": 6, "program": 6, "hope": 6, "help": 6, "fourier": 6, "adder": 6, "multipli": 6, "hhl": 6, "chsh": 6, "model": 6, "blackhol": 6, "evapor": 6}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"quantum": [0, 3, 4, 6], "adder": 0, "multipli": 0, "modular": 0, "sum": 0, "classic": 0, "b": 0, "multipl": 0, "addit": 0, "The": [1, 2], "chsh": [1, 2], "game": 1, "test": 2, "discret": 3, "fourier": 3, "transform": 3, "phase": 3, "estim": 3, "qiskit": 3, "implement": 3, "qpe": 3, "circuit": 4, "model": 4, "blackhol": 4, "evapor": 4, "hhl": 5, "algorithm": 5, "overview": 5, "pre": 5, "process": 5, "hybrid": 5, "topic": 6, "comput": 6}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 56}})