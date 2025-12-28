;;; vaspy-mode.el --- Minor-mode for VASP calculation scripts

;;; Commentary:
;; This mode provides syntax highlighting for VASP INCAR parameters
;; in Python scripts and org-mode files. It highlights common VASP
;; keywords and provides tooltips with brief descriptions.

;;; Code:

(defvar *vasp-keywords*
  '(;; Basic parameters
    "xc" "pp" "encut" "enaug" "ediff" "ediffg"
    "prec" "algo" "nelm" "nelmin"

    ;; K-points and smearing
    "kpts" "gamma" "ismear" "sigma"

    ;; Electronic structure
    "ispin" "magmom" "nupdown" "nbands" "nelect"

    ;; Relaxation and dynamics
    "ibrion" "isif" "nsw" "potim" "ediffg"
    "nfree" "tein" "tebeg" "teend" "smass"

    ;; Output control
    "lwave" "lcharg" "lvtot" "lvhar" "lelf"
    "lorbit" "nedos" "emin" "emax"

    ;; DFT+U
    "ldau" "ldautype" "ldauu" "ldauj" "ldaul"
    "ldauprint" "lmaxmix"

    ;; Van der Waals
    "ivdw" "luse_vdw" "aggac" "param1" "param2"
    "zab_vdw" "bparam"

    ;; Hybrid functionals
    "lhfcalc" "hfscreen" "aexx" "aggax" "aldac"
    "nkred" "time" "precfock"

    ;; GW and BSE
    "algo" "nomega" "nbands" "encutgw"

    ;; Spin-orbit coupling
    "lsorbit" "lnoncollinear" "saxis"

    ;; Density functional
    "gga" "metagga" "lbeefens"

    ;; Dipole corrections
    "ldipol" "idipol" "epsilon" "dipol"

    ;; Symmetry
    "isym" "symprec"

    ;; Performance
    "ncore" "npar" "kpar" "lplane" "lscalu"
    "lreal" "ropt" "addgrid"

    ;; Machine learning
    "ml_lmlff" "ml_istart" "ml_mb" "ml_mconf"

    ;; NEB
    "images" "spring" "lclimb" "ltangentold"

    ;; Common ASE parameters
    "setups" "kspacing" "kgamma")
  "List of common VASP INCAR parameters.")

(defvar *vasp-keyword-docs*
  '(("xc" . "Exchange-correlation functional (PBE, LDA, etc.)")
    ("encut" . "Plane-wave energy cutoff (eV)")
    ("kpts" . "K-point mesh (nx, ny, nz)")
    ("ismear" . "Smearing method (-5=tetrahedron, 0=Gaussian, 1=MP)")
    ("sigma" . "Smearing width (eV)")
    ("ibrion" . "Optimizer (1=quasi-Newton, 2=CG)")
    ("isif" . "What to relax (2=ions, 3=ions+cell)")
    ("nsw" . "Max ionic steps")
    ("ediffg" . "Force convergence (negative = forces in eV/Å)")
    ("ispin" . "Spin polarization (1=non-spin, 2=spin)")
    ("magmom" . "Initial magnetic moments")
    ("ldau" . "Enable DFT+U")
    ("ldautype" . "DFT+U type (1=Liechtenstein, 2=Dudarev)")
    ("ldauu" . "Hubbard U values (eV)")
    ("ivdw" . "Van der Waals correction method")
    ("lhfcalc" . "Enable hybrid functional")
    ("hfscreen" . "Screening parameter for HSE")
    ("lwave" . "Write WAVECAR")
    ("lcharg" . "Write CHGCAR")
    ("nelm" . "Max electronic steps")
    ("ediff" . "Electronic convergence (eV)")
    ("prec" . "Precision (Low, Medium, High, Accurate)")
    ("algo" . "Electronic minimization algorithm"))
  "A-list of brief documentation for VASP keywords.")

(defvar *vasp-keywords-regex*
  (concat "\\(?1:\\<" (regexp-opt *vasp-keywords*) "\\)\\s-*=")
  "Regexp for VASP keywords in Python assignments.")

(defun vasp-tooltip-1 (_win _obj position)
  "Get the tooltip for the keyword under POSITION."
  (save-excursion
    (goto-char position)
    (let* ((word (thing-at-point 'word t))
           (doc (assoc word *vasp-keyword-docs*)))
      (if doc
          (cdr doc)
        "VASP INCAR parameter"))))

(defun next-vasp-keyword (&optional limit)
  "Find VASP keywords up to LIMIT."
  (while (re-search-forward *vasp-keywords-regex* limit t)
    (when (fboundp 'flyspell-delete-region-overlays)
      (flyspell-delete-region-overlays (match-beginning 0)
                                       (match-end 0)))
    (add-text-properties
     (match-beginning 1)
     (match-end 1)
     (list
      'help-echo 'vasp-tooltip-1
      'face 'font-lock-keyword-face
      'mouse-face 'highlight))))

(define-minor-mode vaspy-mode
  "Minor mode for VASP scripts in Emacs.
Provides syntax highlighting for VASP INCAR parameters in Python
scripts and org-mode files."
  :global t
  (if vaspy-mode
      (font-lock-add-keywords
       nil
       '((next-vasp-keyword 1 font-lock-keyword-face)) t)
    (font-lock-remove-keywords
     nil
     '((next-vasp-keyword 1 font-lock-keyword-face))))
  (when (fboundp 'font-lock-flush)
    (font-lock-flush)))

;; Automatically enable for Python and org-mode
(add-hook 'org-mode-hook (lambda () (vaspy-mode 1)))
(add-hook 'python-mode-hook (lambda () (vaspy-mode 1)))

(provide 'vaspy-mode)

;;; vaspy-mode.el ends here
