;;; vaspy-mode.el --- Minor-mode for VASP calculation scripts

;;; Commentary:
;;

;;; Code:

(defun vasp-keywords ()
  "Return a list of vasp keywords."
  (read
   (shell-command-to-string
    "python -c \"from vasp.validate import keywords; print keywords()\"")))

(defvar *vasp-keywords-regex*
  (regexp-opt (vasp-keywords) 'symbols)
  "Regexp for vasp keywords.")

(defun vasp-tooltip-1 (_ _ position)
  "Get the one line tooltip for the keyword under POSITION."
  (save-excursion
    (goto-char position)
    (car (s-split "\n" (vaspdoc (thing-at-point 'word))))))

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
       (match-beginning 0)
       (match-end 0)
       (list
	'help-echo 'vasp-tooltip
	'local-map map
	'face 'font-lock-keyword-face
	'mouse-face 'highlight)))))


(define-minor-mode vaspy-mode
  "Minor mode for Vasp scripts in Emacs."
  :global t
  (font-lock-add-keywords
   nil
   `((next-vasp-keyword 0 font-lock-keyword-face)) t))

(provide 'vaspy-mode)

(provide 'vaspy-mode)

;;; vaspy-mode.el ends here
