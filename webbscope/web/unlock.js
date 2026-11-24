"use strict";
// If the server job dies during a read, conversion or export, the page's buttons stay greyed out. This frees them; the page and its edits are not touched.
(() => {
  const b = document.createElement("button");
  b.textContent = "Job died? Unlock buttons";
  b.title = "Use this if the server job died while a button was greyed out. Your edits on this page are kept. Restart the job first, then press the button again.";
  b.style.cssText = "position:fixed;right:12px;bottom:12px;z-index:10000;font-size:12px";
  b.onclick = () => {
    document.querySelectorAll("#go, #scan, #open, #mergedbox button").forEach((x) => { x.disabled = false; });
    document.querySelectorAll("button[data-d]").forEach((x) => { x.disabled = false; x.textContent = "Draw…"; });
    if (typeof tisState === "function") tisState();             // the tissue-box export stays grey until the boxes are built
    const m = document.querySelector("#tisoutmsg");
    if (m && m.textContent.startsWith("writing")) m.textContent = "";
  };
  document.body.append(b);
})();
