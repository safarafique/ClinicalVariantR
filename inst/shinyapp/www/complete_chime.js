/* ClinicalVariantR — short chime when analysis completes */
(function () {
  function playChime() {
    try {
      var AudioCtx = window.AudioContext || window.webkitAudioContext;
      if (!AudioCtx) return;
      var ctx = new AudioCtx();
      function tone(freq, start, dur) {
        var o = ctx.createOscillator();
        var g = ctx.createGain();
        o.type = "sine";
        o.frequency.value = freq;
        g.gain.setValueAtTime(0.0001, ctx.currentTime + start);
        g.gain.exponentialRampToValueAtTime(0.16, ctx.currentTime + start + 0.02);
        g.gain.exponentialRampToValueAtTime(0.0001, ctx.currentTime + start + dur);
        o.connect(g);
        g.connect(ctx.destination);
        o.start(ctx.currentTime + start);
        o.stop(ctx.currentTime + start + dur + 0.02);
      }
      tone(660, 0, 0.18);
      tone(880, 0.16, 0.28);
    } catch (e) {
      /* autoplay policies / missing Web Audio — ignore */
    }
  }

  function register() {
    if (!window.Shiny || !Shiny.addCustomMessageHandler) return;
    Shiny.addCustomMessageHandler("clinicalvariantrCompleteSound", function () {
      playChime();
    });
  }

  if (window.Shiny) {
    register();
  } else {
    document.addEventListener("shiny:connected", register);
  }
})();
