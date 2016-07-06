;;; vaspy-mode.el --- Minor-mode for VASP calculation scripts

;;; Commentary:
;; The primary purpose of this mode is to provide syntax highlighting
;; for vasp keywords defined in vasp.validate. The syntax highlighting
;; provides tooltips of the first line of documentation for the
;; defined keywords, and makes them clickable to show the whole
;; docstring.

;;; Code:


(defun vasp-keywords ()
  "Return a list of vasp keywords."
  (read
   (shell-command-to-string
    "python -c \"from vasp.validate import keywords; print keywords()\"")))


(defvar *vasp-keywords-regex*
  ;; in hy keywords start with : and in Python they sometimes end in =
  (concat "\\(?1::?\\<" (regexp-opt (vasp-keywords)) "\\)=?")
  "Regexp for vasp keywords.")

(defun vasp-keyword-alist ()
  "Return a list of vasp keywords."
  (read
   (shell-command-to-string
    "python -c \"from vasp.validate import keyword_alist; print keyword_alist()\"")))

(defvar *vasp-keyword-alist*
  (vasp-keywords-alist)
  "A-list of docstrings for vaspy keywords.
Computed on startup.")

(defun vasp-tooltip-1 (_win _obj position)
  "Get the one line tooltip for the keyword under POSITION."
  (save-excursion
    (goto-char position)
    (cadr (assoc (thing-at-point 'word) *vasp-keyword-alist*))))

(defun next-vasp-keyword (&optional limit)
  "Find vasp keywords up to LIMIT."
  (while (re-search-forward *vasp-keywords-regex* limit t)
    (flyspell-delete-region-overlays (match-beginning 0)
				     (match-end 0))
    (let ((map (make-sparse-keymap)))
      (define-key map [mouse-1]
	(lambda ()
	  (interactive)
	  (pydoc (format "vasp.validate.%s" (thing-at-point 'word)))))
      (add-text-properties
       (match-beginning 1)
       (match-end 1)
       (list
	'help-echo 'vasp-tooltip-1
	'local-map map
	'face 'font-lock-keyword-face
	'mouse-face 'highlight)))))


(define-minor-mode vaspy-mode
  "Minor mode for Vasp scripts in Emacs."
  :global t
  (font-lock-add-keywords
   nil
   `((next-vasp-keyword 1 font-lock-keyword-face)) t))


(add-hook 'org-mode-hook (lambda () (vaspy-mode)))
(add-hook 'python-mode-hook (lambda () (vaspy-mode)))

(provide 'vaspy-mode)

;;; vaspy-mode.el ends here
